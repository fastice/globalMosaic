# globalMosaic

Global mosaics of NISAR products, built with the GrIMP `geomosaic` binary: **backscatter**
from GCOV, **coherence** from GUNW, and the **dual-pol ratio** of a like-pol and a cross-pol
backscatter mosaic, on either a lat/lon grid (for Google Earth) or an equal-area projected
grid (EASE-Grid 2.0, Equal Earth, Wagner VII).

Granules are read where they sit at ASF over `/vsicurl/`; nothing is downloaded in bulk.

---

## Install

Needs Python 3.10+, GDAL with its Python bindings, `rasterio`, `shapely`, `pyproj`,
`matplotlib`, `numpy`, `h5py`, `asf_search`, and optionally `cartopy` (coastlines).

```bash
git clone git@github.com:fastice/globalMosaic.git
cd globalMosaic
pip install -e .
installGeomosaic --jobs 8          # clones and builds the six GrIMP C repos
```

`installGeomosaic` fetches `clib`, `cRecipes`, `triangle`, `gdalIO`, `landsatSource64` and
`mosaicSource` into `$HOME/progs/GIT64` and builds `geomosaic` into `$HOME/bin/$(uname -m)`.
It ends with two capability checks, because a binary built from the wrong branch fails only
once a long run is under way. `--gitroot DIR` if the sources live elsewhere; it refuses to
clone a second copy when a populated tree already exists.

**Earthdata credentials** go in `~/.netrc`:

```
machine urs.earthdata.nasa.gov login YOUR_USER password YOUR_PASSWORD
```

The tools set `GDAL_HTTP_NETRC` and a cookie jar themselves; the redirect needs both.
**Without `~/.netrc` every granule read fails and the mosaic comes out empty but valid** —
georeferenced tiles full of no-data, which reads as missing acquisitions rather than as an
authentication problem. `test/checkAccess.sh` checks it before a long run.

Catalogue searches go through `asf_search` (no credentials needed for the search itself);
only reading the granules needs the login.

---

## The three products

### Backscatter (GCOV → gamma0)

```bash
# 1. search a cycle and tile it
globalGCOVTiles --cycle 030 --out cycle30 --direction A --allowV      # lat/lon tiling
easeTiles --catalogue cycle30/catalogue.geojson --hemisphere global \
          --crs ESRI:54076 --tileKm 600 --direction A --allowV \
          --out run/tiles                                             # projected tiling

# 2. mosaic every tile (one geomosaic process per job)
runGeomosaicTiles run/tiles --work run --granules vsicurl --psResM 1000 --nProc 32

# 3. look at it
globeFigure run/tiles/gamma0 --out globe.png --width 6000 --coast
```

Tiles are int16 **dB x 100**, no data −3000.

### Coherence (GUNW → coherenceMagnitude)

```bash
easeTiles --catalogue cycle24gunw/catalogue.geojson --hemisphere global \
          --crs ESRI:54076 --tileKm 600 --direction A --allowV --allGranules \
          --out run/tiles
gunwCoherence run/tiles --work run --psResM 1000 --nProc 32
globeFigure run/tiles/coherence --out coherence.png --width 6000 --coast
```

Tiles are Float32 **0–1 linear**, no data NaN. `gunwCoherence` writes a ~1.5 kB VRT per
granule that points at the coherence array at ASF and lets `geomosaic -geo` read windows from
it; nothing is copied locally.

### Dual-pol ratio (GCOV like-pol / cross-pol)

openSEPPO's `-dpratio` recipe: like-pol over cross-pol, `DH` -> HHHH/HVHV, `DV` -> VVVV/VHVH,
`QP` -> HHHH/HVHV.

```bash
bash test/runCycleRatio.sh 033                      # both bands, the ratio and the figure
MODE=amp bash test/runCycleRatio.sh 033
```

Because geomosaic writes gamma0 as int16 dB x 100, the ratio is arithmetic on two dB mosaics,
so `ratioMosaic` needs nothing from the C:

| `--mode` | value | type | no data |
|---|---|---|---|
| `dB` (default) | `likeDb - crossDb` | Float32 | NaN |
| `pwr` | `10 ** (d / 10)` | Float32 | NaN |
| `amp` | `10 ** (d / 20) * 1000`, clipped [1, 65535] | UInt16 | 0 |
| `DN` | `DN(like) - DN(cross) + 127`, clipped [1, 255] | Byte | 0 |
| `rgb` | R = DN(like), G = DN(cross), B = the DN ratio | 3 x Byte | 0 |

A pixel is no data where **either** input is, and for `amp`/`DN` that is decided *before*
clipping, so a saturated pixel saturates instead of vanishing. DN follows openSEPPO's
`power_to_dn_uint8` -- `dB = -31.15 + DN * 0.15`, clamped, DN 0 reserved -- and `DN`/`rgb`
difference the CLAMPED DNs, which is not the same as scaling the dB difference: cross-pol
below -31 dB pins to DN 1, routine over desert and water.

`rgb` is the published false-colour composite, and for it the **like-pol input should be the
FULL-cover mosaic**: where there is no cross-pol, all three bands take DN(like) so single-pol
ground -- the ice sheets, and the 16% of the globe acquired 5 MHz VV -- renders as greyscale
backscatter instead of dropping out.

`--minCross DB` drops pixels whose cross-pol is at or below a threshold. Once the denominator
is at the noise floor (NESZ is about -25 dB for NISAR L-band) the ratio is a lower bound on a
quantity that was not measured. Open water is the case that shows it: smooth water barely
depolarises, HV sits at the floor, and the ratio comes out as the brightest thing in the
scene. `-28` is a reasonable value; in `rgb` those pixels fall back to grey, like single-pol
ground, so backscatter is still shown where the ratio is not.

`--stretch LO HI` applies a per-band percentile stretch (off by default, which writes raw DN).
Single-pol pixels take the RED band's limits in all three bands, or the greyscale fallback
would come out tinted. There is no published stretch to match: the reference map is scaled by
**SEPPO**, a separate package from openSEPPO, and its scaling is not documented. Note also
that cross-pol is always weaker than like-pol, so `B` is confined to [127, 255] and blue never
falls below mid scale -- a raw DN composite cannot render forest as pale yellow.

**The two bands must come from the same granule at every pixel**, or the ratio compares two
acquisitions. That is what `easeTiles --ratioCover` (like-pol) and `--crossPol` (cross-pol)
are for: both restrict the cover to granules that *have* a cross-pol term, dropping `SH` and
`SV`, and the cover is deterministic given the same candidates, so the two runs agree. Use
`--allowV` as well or `DV` granules never reach the cover.

**An existing like-pol mosaic cannot be reused for a ratio.** Built without `--ratioCover`,
its cover includes single-pol granules; the cross-pol pass skips those and fills the pixels
from a different acquisition, and the result looks fine exactly where it is wrong. Expect the
ratio to cover less ground than a plain backscatter mosaic too -- about 72% of cycle-30
granules carry a cross-pol term.

---

## Ready-made drivers (`test/`)

| script | what it builds |
|---|---|
| `runEaseMosaic.sh` | one global backscatter mosaic from every cycle (fill catalogues), projected |
| `runCycleGCOV.sh <cycle>` | backscatter for ONE cycle, minimal cover |
| `runCycleCoherence.sh <cycle>` | coherence for ONE cycle, every overlapping pair averaged |
| `runCycleRatio.sh <cycle>` | dual-pol ratio for ONE cycle: like-pol + cross-pol on one cover, then the ratio |
| `rebuildAllGCOV.sh` / `rebuildAllCoherence.sh` | wipe and rebuild a series of cycles |
| `makeGoogleEarth.sh`, `runFull.sh`, `runSmoke.sh` | the lat/lon Google Earth product |

Smaller helpers in the same directory: `setup.sh` and `checkAccess.sh` (environment and ASF
credentials), `monitor.sh` and `netMonitor.sh` (progress and download rate), `gaps.sh`,
`histograms.sh`, `runStatistics.sh`, `statsTest.sh` (analysis of a finished run),
`runCoherence.sh` (the lat/lon coherence path), `runAfter.sh` (chain a job behind another),
and `restartGoogleEarth.sh`, `retitleGoogleEarth.sh`, `tarGoogleEarth.sh` (maintain and ship
the Google Earth product).

Each takes its settings from the environment: `ROOT`, `CRS`, `TILEKM`, `RES`, `NPROC`,
`DIRECTION`, `CYCLES`. A per-cycle series lands on one grid — the tiling depends only on
`CRS` and `TILEKM`, never on the data — so cycles difference directly, and so do backscatter
and coherence of the same cycle.

---

## Projections

The output grid is whatever `geomosaic` can be given, so any CRS PROJ can parse works —
EPSG code, ESRI code or proj string. All of these were run against the binary:

| kind | example | how it is handled |
|---|---|---|
| geographic | `EPSG:4326` | grid units are DEGREES, not metres |
| polar stereographic | `EPSG:3031`, `3413`, `3995`, custom `+proj=stere` | the legacy exact path; conformal |
| UTM | `EPSG:32610`, `32633`, `32760` | through PROJ; conformal |
| equal area with a code | `EPSG:6933` (EASE Global), `6931`/`6932` (EASE polar), `3338`, `5070` (Albers) | positions only — see below |
| equal area without a code | Equal Earth, `ESRI:102001`, any proj string | same, passed with `-wkt` instead of `-epsg` |
| Wagner VII | `ESRI:54076` | special: PROJ has no inverse, so `grimpProj.c` carries its own |

**"Positions only"** is what `geomosaic` prints for a non-conformal CRS, and for a
backscatter or coherence mosaic it is simply a note: a scalar is resampled, so an equal-area
or any other projected grid is fine. It matters for **velocity**: `mosaic3d` decomposes into
vx/vy with one grid angle and one isotropic scale, which only holds on a conformal grid, so
it refuses anything but polar stereographic or UTM.

Two hard limits, both of which error out rather than producing a wrong answer:

- **WGS84 only.** A CRS on another ellipsoid is refused — `EPSG:3411` (Hughes 1980,
  semi-major 6,378,273) fails with the axis quoted. GRS80 passes, being WGS84's axis to
  5e-9.
- **No false easting or northing on a polar stereographic**; the exact legacy path cannot
  express one.

In this package, `easeTiles --crs` takes any of these for `--hemisphere global`: a CRS with
an EPSG code reaches `geomosaic` as `-epsg`, one without (Equal Earth, Wagner VII) is written
to `proj.wkt` and passed as `-wkt`. The `-wkt` path and Wagner VII need a `geomosaic` built
from current `mosaicSource`; `-epsg` works with any build.

An ESRI code must carry its prefix (`ESRI:54076`) — `-epsg 102001` fails, because that is
not an EPSG code.

---

## Every command

| command | what it does |
|---|---|
| `globalGCOVTiles` | search ASF for a cycle and tile it in lat/lon (6x6 deg, polar caps beyond 84 deg); tiered greedy cover 77 → 40 → 20 → 5 MHz |
| `easeTiles` | the same for a projected grid: EASE-Grid 2.0 polar (6931/6932), or `--hemisphere global` with `--crs` (Equal Earth, EPSG:6933, Wagner VII) |
| `scanMasks` | valid-DATA footprints of GCOV (from the mask) or GUNW (where coherence is finite), for `--validFootprints` |
| `runGeomosaicTiles` | run `geomosaic` per tile, crop the feather margin, build the VRT hierarchy, draw a live progress picture |
| `gunwCoherence` | coherence tiles via `geomosaic -geo`, on any tiling including ones PROJ cannot invert |
| `ratioMosaic` | dual-pol ratio of a like-pol and a cross-pol mosaic (`dB`/`pwr`/`amp`/`DN`), streamed blockwise |
| `mosaicCoherence` | the older lat/lon coherence path (warps each granule; also builds Google Earth) |
| `frameStatistics`, `mosaicStatistics` | temporal mean/sigma/CV per frame, then mosaicked with `-geo` |
| `globeFigure` | the mosaic as a globe on filled ocean, in its own projection |
| `mosaicGaps`, `gapReport`, `gapsFromPicture` | find holes: against granule footprints (projected grids) or a land mask (lat/lon) |
| `coverageMaps`, `makeQuickLook`, `backscatterHistograms`, `monitorRun` | coverage, quick looks, statistics, progress |
| `makeGoogleEarth`, `makeGoogleEarthTiered`, `oceanAltitude`, `linkOnly` | KML superoverlays |
| `installGeomosaic`, `benchmarkASFRead` | build the C binary; measure ASF read throughput |
| `notAcquired` | ground no granule covers |

---

## Things that will bite you

- **Direction, bandwidth and polarization come from the file name**, not from ASF metadata:
  `flightDirection` is wrong for about 4% of granules.
- **5 MHz is almost all V-only.** Ascending 5 MHz granules across cycles 023–031: 2,472
  without `--allowV`, 20,771 with it. The 18,299 you lose are most of the Sahara and the sea
  ice.
- **An acquisition footprint overstates.** A frame's valid data can stop ~100 km short of its
  footprint. `scanMasks` + `--validFootprints` makes the cover choose on real data.
- **A footprint that crosses 180 must be split before projecting** — `easeTiles` does this. It
  is not optional: an unsplit one spans the whole map and makes the cover believe every tile
  at that latitude is already covered, leaving bands right across the product.
- **Wagner VII has no inverse in PROJ.** `geomosaic` carries its own (`grimpProj.c`), and the
  Python side has a matching `wagner7Inverse`. `gdal.Warp` into such a CRS fails; use
  `gdal.Translate` for a plain resize.
- **`CPL_VSIL_CURL_CHUNK_SIZE` and the vsicurl CACHE decide how much you re-download.** A GCOV
  band is gzipped in 512x512 chunks (0.26 MB compressed each), and geomosaic walks it in
  strips, so a large curl chunk fetches far more than the chunk it wanted. One HHHH band is
  1.60 GB stored; reading it for a global tile over `/vsicurl`:

  | chunk | cache | transferred | time |
  |---|---|---|---|
  | 4 MB | 16 MB (GDAL default) | 12760 MB | 1055 s |
  | 16 kB (GDAL default) | 16 MB | 3475 MB | 699 s |
  | 64 kB | 256 MB | 1795 MB | 317 s |

  The cache is the half that is easy to miss: one chunk ROW is 20.8 MB against a 16.4 MB
  default cache, and a strip straddles two chunk rows, so each row is evicted before its
  second use and downloaded twice. `runGeomosaicTiles` now defaults to 64 kB / 256 MB (per
  process, so 32 jobs hold 8 GB) and both are overridable from the environment. A 4 MB chunk
  remains right for one large sequential stream, which is what `benchmarkASFRead.py` measures
  and is not how a GCOV is read.

- **`GDAL_CACHEMAX` defaults to 5% of PHYSICAL RAM, per process** -- 6.2 GB on a 124 GB
  machine, 25 GB on a 503 GB one -- and you run one process per `--nProc`, so the default is
  an out-of-memory waiting for a big job count. Both drivers now set 256 MB, overridable:
  enough for the ~158 MB of blocks one strip read touches, and 20 GB across 80 jobs. geomosaic
  reads each strip once and never revisits it, so a bigger cache only holds blocks nothing
  will ask for again.

- **The NISAR GDAL plugin trades bytes for latency.** With `gdal_NISAR.so` installed,
  `GRIMP_NISAR_DRIVER=remote` sends the big 2-D arrays through it for `/vsi` paths: the same
  band reads in ~100 s instead of 317 s, bit-identical, because it keeps ~30 range requests in
  flight while libhdf5 blocks on one at a time (a flat ~0.2 s each). It moves slightly MORE
  data (1.24x the stored bytes against 1.10x), so the win is latency, not volume -- it pays
  when bandwidth is idle and stops paying once `--nProc` jobs saturate the link, where the
  extra bytes and ~30 sockets per job turn into a liability. Off by default; use it for one or
  a few jobs, not a 32-way global run. It is 6-10x SLOWER on local files and never used for
  them.
- **A ratio needs its own cover.** Both bands must come from one granule per pixel, so the
  like-pol half of a ratio pair is built with `easeTiles --ratioCover`, not reused from a
  normal backscatter run. Single-pol granules (`SH`, `SV`) have no cross-pol term at all.
- **Per-cycle series take one cycle only.** `--fillCatalogues` closes gaps by borrowing other
  cycles, which is right for a single best-coverage product and wrong for a time series.

## Tests

`test/runSmoke.sh` is the quickest end-to-end check. `test/checkAccess.sh` verifies ASF
credentials and `geomosaic` before a long run.
