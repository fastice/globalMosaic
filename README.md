# globalMosaic

Global mosaics of NISAR products, built with the GrIMP `geomosaic` binary: **backscatter**
from GCOV and **coherence** from GUNW, on either a lat/lon grid (for Google Earth) or an
equal-area projected grid (EASE-Grid 2.0, Equal Earth, Wagner VII).

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

## The two products

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

---

## Ready-made drivers (`test/`)

| script | what it builds |
|---|---|
| `runEaseMosaic.sh` | one global backscatter mosaic from every cycle (fill catalogues), projected |
| `runCycleGCOV.sh <cycle>` | backscatter for ONE cycle, minimal cover |
| `runCycleCoherence.sh <cycle>` | coherence for ONE cycle, every overlapping pair averaged |
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
- **`CPL_VSIL_CURL_CHUNK_SIZE` cuts both ways**: 4 MB is best for geomosaic's large sequential
  reads, and 44x worse than the default for the scattered small reads that open an HDF5.
- **Per-cycle series take one cycle only.** `--fillCatalogues` closes gaps by borrowing other
  cycles, which is right for a single best-coverage product and wrong for a time series.

## Tests

`test/runSmoke.sh` is the quickest end-to-end check. `test/checkAccess.sh` verifies ASF
credentials and `geomosaic` before a long run.
