# globalMosaic

Global, lat/lon-tiled NISAR GCOV backscatter mosaics built with GrIMP `geomosaic`.

| command | what it does |
|---|---|
| `globalGCOVTiles` | search ASF for one cycle (routine `L2_PR` products only), tile the area in lat/lon (6 x 6 deg to 60 deg, wider bands toward the poles, one polar stereographic tile per cap beyond 84 deg), and pick a non-redundant covering set of granules per tile (tiered 77 -> 40 -> 20 -> 5 MHz greedy cover); writes `tiles.geojson`, per-tile CSVs and geomosaic `-gcov` yamls |
| `coverageMaps` | coverage map of a tiling run: mid-latitudes plus Arctic and Antarctic polar views, by bandwidth, uncovered land in red |
| `makeGoogleEarth` | EPSG:4326 KML superoverlay of compressed 8-bit PNG tiles for Google Earth (gdal2tiles geodetic, `--zoom 0-9` default = ~305 m pixels; 11 = full 3" resolution), polar caps warped in |
| `makeQuickLook` | quick-look PNGs of a run: one per tile (0.01 deg), a global one (0.05 deg) and one per polar cap, in `<work>/quicklooks/` |
| `makeGoogleEarthTiered` | Google Earth with the zoom limit following the data: base 0-9 everywhere (ocean and sea ice stop there), land 5-10 (2x averaged + 3x3 smoothed, no-data aware), detail 7-11 over 40/77 MHz land, ice-sheet margins (200 km of the grounded coast, plus ice shelves) and glaciers; one doc.kml loads all three |
| `runGeomosaicTiles` | run `geomosaic` on every tile (one thread per process, many processes; `-calOutput gamma0`, `dem none`, `-epsg 4326` / `3031` / `3413`), crop the feather margins and build the VRT hierarchy tile -> latitude band -> global |

Direction, bandwidth and polarization come from granule file names (ASF's `flightDirection`
is wrong for about 4% of granules). Needs a geomosaic with geographic output and `dem none`.

`cycle30/` holds the cycle-030 tile descriptions: `catalogue.geojson` (the ASF search) and the
tilings for `both` directions, `ascending` and `descending` (`--allowV --exclude5MHzOcean`).

## Running (test/)

Run files are yamls whose keys are `runGeomosaicTiles` option names; the command line overrides
them and the resolved settings are saved to `<work>/run.yaml`.

```bash
test/setup.sh            # pip install, build geomosaic, check ~/.netrc
test/runSmoke.sh [test/smokeAsc.yaml]   # default test/smokeAsc.yaml (ascending); test/smoke.yaml = both; work: says where: 2 tiles + Google Earth, ~10 min
test/runFull.sh [test/fullDesc.yaml]    # default test/fullAsc.yaml (ascending); test/full.yaml = both: 930 jobs, 24 procs, then Google Earth (background)
test/monitor.sh [yaml]   # progress picture from disk, safe while a run goes: <work>/quicklooks/monitor.png
test/makeGoogleEarth.sh  # redo the Google Earth tiles (other stretch/zoom)
```

## Pixel convention

GrIMP C programs (geomosaic included) work in **pixel-centre** coordinates: the input-file header
origin is the centre of the first pixel. GeoTIFFs and VRTs are **edge** based; geomosaic's writer
converts (origin minus half a pixel) unconditionally. `runGeomosaicTiles` therefore writes headers
as tile edge + half a pixel, and crops using each output's own (edge-based) geotransform.
geomosaic's `-center` flag is a no-op; the convention does not depend on it.
