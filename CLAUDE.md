# CLAUDE.md — globalMosaic

Guidance for an assistant asked to build a NISAR mosaic with this package. Read the README
for what each command does; this file is how to turn a request into the right commands.

A request is usually "a mosaic of **what** for **when**, over **where**, in **which
projection**". Resolve those four, then follow the matching recipe.

## 0. Access — check this first

Granules are never downloaded in bulk: GDAL reads windows straight from ASF over
`/vsicurl/`. That needs an Earthdata login in **`~/.netrc`**:

```
machine urs.earthdata.nasa.gov login YOUR_USER password YOUR_PASSWORD
```

Without it every read fails and a run produces correctly georeferenced, entirely empty tiles
— which looks like missing data, not like an auth problem. Check it before a long run:

```bash
bash test/checkAccess.sh
```

Searching is separate: catalogues come from **`asf_search`** (`globalGCOVTiles
--cycle`, with `CMR_TIMEOUT` raised to 600 s because a whole-cycle query outgrows the
30 s default). The search itself needs no credentials; reading the data does. The `s3`
granule mode additionally fetches temporary S3 credentials from ASF using the same
`~/.netrc`, and those expire after an hour, so `vsicurl` is the safer default for long runs.

---

## 1. What — backscatter or coherence

| they say | product | source | tiles |
|---|---|---|---|
| backscatter, gamma0, sigma0, "an image", "how it looks" | **GCOV** | `runGeomosaicTiles` | int16, dB x 100, no data −3000 |
| coherence, correlation, "how well it correlates", InSAR quality | **GUNW** | `gunwCoherence` | Float32, 0–1 linear, no data NaN |

They differ in one more way that matters. Backscatter takes the **minimal cover** — one good
observation per place, best bandwidth first. Coherence takes `--allGranules`, so every
overlapping pair is averaged and a cell is the mean coherence rather than whichever pair a
cover happened to pick. Keep that distinction unless told otherwise.

## 2. When — a period becomes cycles

NISAR repeats every **12 days**. Cycle N's track 1 is acquired on

```
2026-06-02 + 12 * (N - 22)        # cycle 22 = 2026-06-02, verified against the archive
```

and the cycle runs 12 days from there. Cycles 24–33:

| cycle | acquisitions | | cycle | acquisitions |
|---|---|---|---|---|
| 24 | 2026-06-26 → 07-07 | | 29 | 2026-08-25 → 09-05 |
| 25 | 2026-07-08 → 07-19 | | 30 | 2026-09-06 → 09-17 |
| 26 | 2026-07-20 → 07-31 | | 31 | 2026-09-18 → 09-29 |
| 27 | 2026-08-01 → 08-12 | | 32 | 2026-09-30 → 10-11 |
| 28 | 2026-08-13 → 08-24 | | 33 | 2026-10-12 → 10-23 |

**Coherence needs the pair**, so a cycle is only usable once the *following* cycle's
acquisitions are complete — cycle N is ready 24 days after it starts. Backscatter needs only
the cycle itself. Check `date` before promising the latest cycle.

A period spanning several cycles is a **choice to put to the user**:

- *a series* — one map per cycle, the same grid each time, so they difference. This is what
  you want for change over time. Use `test/runCycleGCOV.sh` or `test/runCycleCoherence.sh`
  per cycle.
- *one map* — best coverage over the whole period, with `--fillCatalogues` filling each
  tile's gaps from other cycles. Good for a single complete picture, wrong for a time series
  because each place comes from a different date.

## 3. Where — an area becomes a tile list

The tiling always covers the globe; a regional mosaic is the same tiling with only some tiles
built. Do not build a special grid for an area — the grid depends only on `--crs` and
`--tileKm`, which is what lets any two runs be compared pixel for pixel.

```bash
tilesFor run/tiles --bbox -75 -60 -30 15            # what covers South America
runGeomosaicTiles run/tiles --work run --granules vsicurl --psResM 1000 --nProc 32 \
    --tiles $(tilesFor run/tiles --bbox -75 -60 -30 15 --names)
```

`--point LON LAT` gives the single tile holding a place. `gunwCoherence` takes `--tiles` too.

## 4. Which projection

| ask | `--crs` | note |
|---|---|---|
| "equal area", "global", default | `ESRI:54076` (Wagner VII) | lens-shaped whole Earth; **needs a geomosaic built from current mosaicSource** |
| cylindrical, "like a normal map" | `EPSG:6933` | EASE-Grid 2.0 Global; works with any geomosaic |
| "Equal Earth" | `+proj=eqearth +lon_0=0 +datum=WGS84 +units=m +no_defs` | needs the same current geomosaic |
| one pole | `--hemisphere north` / `south` | EASE-Grid 2.0 polar, EPSG 6931/6932 |
| Google Earth, KML | — | different path: `globalGCOVTiles` + `makeGoogleEarthTiered` |

Wagner VII and Equal Earth have no EPSG code and reach geomosaic as `-wkt`; Wagner VII also
has no PROJ inverse, which only the current `mosaicSource` handles. If `installGeomosaic` has
not been run on this machine, prefer `EPSG:6933` or run it first.

---

## Recipes

**One cycle, global, backscatter** (the common request):

```bash
NPROC=32 bash test/runCycleGCOV.sh 030
```

**One cycle, one region:**

```bash
easeTiles --catalogue cycle30/catalogue.geojson --hemisphere global --crs ESRI:54076 \
          --tileKm 600 --direction A --allowV --out run/tiles
runGeomosaicTiles run/tiles --work run --granules vsicurl --psResM 1000 --nProc 32 \
    --tiles $(tilesFor run/tiles --bbox 10 35 45 72 --names)
globeFigure run/tiles/gamma0 --out europe.png --width 4000 --coast
```

**A coherence series over several cycles:**

```bash
for c in 028 029 030; do bash test/runCycleCoherence.sh $c; done
```

**One best-coverage map over a long period:** `test/runEaseMosaic.sh`, which passes
`--fillCatalogues` for every cycle it has.

A catalogue it does not have is searched from ASF automatically (~10 min). Catalogues for
cycles 023–031 ship in `cycle30/`, GUNW ones in `cycle30gunw/` and `cycle24gunw/`.

---

## Rules that are not optional

- **`--allowV`.** 5 MHz is almost entirely V-polarization: 2,472 ascending 5 MHz granules
  without it, 20,771 with. Leaving it off empties the Sahara and the sea ice.
- **One direction per product.** Ascending and descending are separate mosaics; never mix
  them in one, and give a descending run its own `ROOT` or it overwrites the ascending one.
- **One run per work directory.** Two concurrent runs sharing `--work` corrupt tiles. Stop
  one first; a restart is cheap because finished tiles are skipped.
- **`rebuildAll*.sh` deletes.** To *add* a cycle to an existing series call
  `runCycleGCOV.sh` / `runCycleCoherence.sh` directly.
- **Changing `--crs` or `--tileKm` means a new tiling directory.** Reusing one across
  projections silently produces nonsense.
- **Re-tiling after changing selection flags needs the tiles rebuilt.** `runGeomosaicTiles`
  and `gunwCoherence` notice when a tile's granule list changed and redo it; deleting
  `tiles.geojson` is what forces the tiling itself to be recomputed.

## Cost, so you can warn them

At 1 km with 32 concurrent jobs, a global cycle is roughly 900–1,400 tiles at a few minutes
each: **4–8 hours**, dominated by reading granules over the network, not by CPU. A region is
proportionally less. Finer posting scales with pixels — 500 m is 4x the work of 1 km. Always
start long runs under `nohup` and tell the user where the log and `quicklooks/progress.png`
are.

## When a mosaic has holes

`mosaicGaps` reports them against the granule footprints and names the tile each falls in.
Before blaming the data, check: was `--allowV` set; does the tiling predate the current code;
is it a per-cycle product where that ground simply was not acquired. The README's "Things
that will bite you" lists the failure modes that produced real holes here.
