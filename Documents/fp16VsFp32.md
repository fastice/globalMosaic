# float16 vs float32 for NISAR products

Measured 2026-09-28 on helheim, files on `/Volumes/insar1` (local XFS, not NFS), against a real
NISAR RSLC (`NISAR_L1_UR_RSLC_031_126_A_030_2005_DHDH_...`, 4.0 GB, downloaded from ASF).
Reproduce with `bench/fp16bench.py`; build the C readers with `bench/buildFp16read.sh`.

All reads evict the file from the page cache first (`posix_fadvise(POSIX_FADV_DONTNEED)`, no root
needed). Without that, every read after the first measures the page cache and the comparison is
meaningless — a warm-cache row is included to show how large that effect is.

## The finding that reframes everything: `complex64_zero_mantissa`

The RSLC's own ISCE3 runconfig (`metadata/processingInformation/parameters/
runConfigurationContents`) says it outright:

```
data_type: complex64_zero_mantissa
compression_enabled: true
compression_type: gzip
compression_level: 4
shuffle: true
```

This is a deliberate encoding, not an oversight. The low mantissa bits of each float32 are
zeroed so that shuffle+gzip can strip them. Measured: the low **13** bits are zero on all
33,554,432 components, leaving 23 - 13 = **10 mantissa bits -- exactly float16's mantissa**.

It buys three things simultaneously:

| | |
|---|---|
| float16 mantissa precision | 10 bits; the information content really is fp16 |
| **float32 exponent range** | 8 bits, +/-38 decades, against float16's 5 bits and +/-5 |
| universal readability | it is ordinary complex64; every tool reads it |

plus compressibility close to true fp16, because the dead bits are constant.

**This is the dynamic-range argument for their choice over complex32.** This scene's I/Q
components span 189.9 dB in power against float16's 180.6 dB normal range, and 0.0130 % of them
fall below float16's smallest normal value. Under `complex64_zero_mantissa` those samples keep
full relative precision, because the float32 exponent is untouched. Convert to real complex32
and they degrade into subnormals. The penalty is small here -- they are noise-floor samples --
but it is the thing the encoding was designed to avoid.

Consequence for a storage decision: converting c64-zero-mantissa -> c32 buys **1.09x on disk**
and costs the exponent range plus GDAL readability. The saving is already banked in the
encoding; there is very little left to take.

## Corroborating measurements


**NISAR RSLCs are already float16 data stored in float32 containers.**

| | |
|---|---|
| I/Q components exactly representable in float16 | **99.9915 %** |
| low 13 mantissa bits of the float32 that are non-zero | **0 of 33,554,432** (0.0000 %) |
| distinct component values present | **35,064** (float16 has 65,536 codes) |
| the 0.0085 % that are not exact | all below 6.1e-5, i.e. fp16 subnormals, at the noise floor |

Converting complex64 → complex32 on this product gives **median and p99 phase error of exactly
zero**, max 0.0057 deg, and coherence 1.00000000. There is no information to lose, because the
processor already quantised to 16-bit precision before writing 32-bit containers.

Checked on routine (`PR`) as well as urgent (`UR`) products, 20 MHz and 40 MHz, dual-pol: all are
`{r:float32, i:float32}`. **ASF does not ship complex32.** So fp16 for RSLC is a conversion you
would make, not a format you can read natively.

## Bits per sample — and why compression already does fp16's job

| | bits per complex sample |
|---|---|
| HDF5 container `{r:f32, i:f32}` | 64 |
| log2(35,064 distinct values) | 30.2 |
| Shannon entropy of the components | 27.1 |
| **on disk, c64 + gzip-4 + shuffle** | **30.1** |
| **on disk, c32 + gzip-4 + shuffle** | **27.6** |

This is the single most important number here. The shuffle filter groups the always-zero low
mantissa bytes together and gzip removes them, so **c64 already stores ~30 bits per sample, not
64**. Dropping to c32 saves only 2.5 bits — a **1.09x** file-size reduction, not 2x.

The 30.1 bits measured on disk against 30.2 bits of log2(distinct) says gzip+shuffle is already
at the redundancy limit for this data. There is almost nothing left for fp16 to take out.

Measured file sizes (8192 x 8192 window, gzip-4 + shuffle):

| | file | vs fp32 |
|---|---|---|
| float32 (magnitude) | 216.0 MB | — |
| float16 | 105.1 MB | **2.06x smaller** |
| complex64 | 252.7 MB | — |
| complex32 | 231.7 MB | **1.09x smaller** |

Note the difference: the *real* float32 band halves properly, while *complex* barely moves. The
magnitude array was rescaled and is smooth; the raw complex is noise-like in phase, so gzip has
nothing to work with on the fp16 version while it still strips the dead mantissa bytes from the
fp32 one.

**So the storage argument for fp16 depends entirely on whether you compress.** Uncompressed, it
is exactly 2x. With gzip+shuffle it is 2.06x for a smooth real band and 1.09x for raw complex.

## Read speed: fp16 is faster ONLY if you avoid HDF5's converter

8192 x 8192 window, cold cache, seconds:

| reader | float32 | float16 | fp16 vs fp32 |
|---|---|---|---|
| C, GDAL 3.8.4 | 1.29 | 6.15 | **0.21x** (4.8x slower) |
| C, GDAL 3.9.0 | 1.11 | 6.16 | **0.18x** (5.6x slower) |
| C, HDF5 converts (1.10.10) | 1.20 | 6.76 | 0.18x |
| C, HDF5 converts (1.14.3) | 1.17 | 5.55 | 0.21x |
| **C, HDF5 + lookup table (1.10.10)** | 1.43 | **0.89** | **1.61x faster** |
| **C, HDF5 + lookup table (1.14.3)** | 1.37 | **0.88** | **1.56x faster** |
| python GDAL 3.9 | 3.25 | 6.78 | 0.48x |
| h5py, native dtype | 1.08 | 0.53 | 2.03x faster |
| h5py, asking for float32 | 1.30 | 6.91 | 0.19x |
| python GDAL, WARM cache | 3.12 | 7.32 | 0.43x |

**The newer GDAL does not fix this.** 3.8.4 and 3.9.0 are within noise of each other; both are
~5x slower on fp16 than fp32. Anything that lets libhdf5 do the conversion pays the same
penalty, in C or in Python.

Note the warm-cache row: it is no faster than cold. These reads are bound by decompression and
conversion, not by disk — which is why the conversion path dominates everything.

## Why a lookup table, and what it is not

`geoMosaic/gcovMosaic.c` reads the raw 16-bit values using the file's own datatype as the memory
type — so HDF5 copies bits instead of converting — and expands them through a 65,536-entry
table (256 KB, fits in L2).

Pure conversion, in memory, 67 Mpx, no I/O:

| | Mpx/s |
|---|---|
| numpy `astype` (hardware F16C) | 301.3 |
| 65,536-entry lookup table | 242.8 |
| **HDF5's converter, on read** | **10.6** |

So the LUT is **1.24x slower than hardware F16C** — it is not the fastest way to convert. It is
chosen because it needs no architecture-specific compiler flag and still builds on arm64. The
win is not over F16C; it is that **both are ~25x faster than what libhdf5 does internally**. The
bottleneck was never the arithmetic.

## Write speed

gzip-4 + shuffle, 8192 x 8192, h5py:

| | seconds | Mpx/s |
|---|---|---|
| float32 | 7.21 | 9.3 |
| float16 | 3.58 | 18.8 |
| complex64 | 13.55 | 4.9 |
| complex32 | 9.96 | 6.7 |

fp16 writes **2.0x faster** for the real band — writing is dominated by compressing bytes, and
there are half as many. Complex gains less (1.36x) for the same reason its file barely shrinks.

## Precision

On the real data, float32 -> float16 -> float32:

| | |
|---|---|
| median relative error | 1.7e-4 (**0.0007 dB**) |
| p99 | 4.3e-4 (0.0019 dB) |
| max | 6.5e-3 (0.0281 dB) |
| flushed to zero | 0 of 67,108,864 |
| overflowed to inf | 0 |
| distinct levels | 13,320,243 -> 12,724 (**1047x fewer**) |

0.0007 dB against a 0.01 dB output quantisation is irrelevant. The risk is not error magnitude,
it is **sample loss**: values below 5.96e-8 flush to zero and vanish, and values above 65,504
become inf. Neither happened here, but both are silent — no tolerance on the values would catch
them. Gate on `min(valid) > 1e-7` and `max < 65504` per product, not on a dB difference.

## Dynamic range, specifically for SLC

| | |
|---|---|
| float16 normal range | 6.104e-5 .. 65,504 = 1.07e9 = **180.6 dB in power** |
| this scene's I/Q components | 6.6e-9 .. 20.7 = 3.12e9 = **189.9 dB** |
| components that overflow | **0** |
| components in the subnormal region | 4,352 (**0.0130 %**) |
| components that flush to zero | 8 (0.00002 %) |

The scene is marginally wider than float16's normal range, and the excess sits entirely at the
bottom — noise-floor samples, where reduced precision is harmless. Since the data is already
fp16-valued anyway (above), **there is no dynamic-range argument against complex32 for NISAR
RSLC as currently produced**. That could change if a future processor stopped pre-quantising.

## Summary

**Advantages of fp16**
- Halves uncompressed size and in-memory footprint, exactly.
- ~2x faster to write (half as many bytes to compress).
- ~1.6-2x faster to read **if** you bypass HDF5's converter (LUT in C, native dtype in h5py).
- For NISAR RSLC, lossless in practice: the data is already fp16-valued.

**Disadvantages**
- **~5x slower to read through GDAL**, on both 3.8.4 and 3.9.0, and through any path that lets
  libhdf5 convert. This is the big one, and a newer GDAL does not fix it.
- **Compression eats most of the storage win**: 1.09x for raw complex, 2.06x for a smooth real
  band. If you gzip+shuffle anyway, fp32 already stores ~30 bits per sample.
- Silent failure modes at the ends of the range (flush-to-zero, overflow-to-inf) that no value
  tolerance detects.
- 1047x fewer distinct levels; fine for backscatter, think again for anything differenced.
- Not natively readable as complex32 by tools expecting complex64, and `H5T_IEEE_F16LE` only
  exists from HDF5 1.14.4 — neither 1.10.10 nor 1.14.3 has it, so the half type must be built by
  hand (`H5Tset_fields` / `H5Tset_ebias`), which is what the C benchmark and the slim products do.

**When it is worth it**: uncompressed or lightly compressed storage, memory-bound pipelines, and
readers you control (so you can use the LUT). **When it is not**: gzip+shuffle archives read
through GDAL by other people's tools — there you pay a 5x read penalty for a 9 % space saving.
