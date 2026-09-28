# cycle30 tilings

| file | what |
|---|---|
| `catalogue.geojson` | cycle 030 GCOV granules (ASF search), `L2_PR`, one per scene |
| `catalogue029.geojson`, `catalogue031.geojson` | cycles 029 and 031, used only to fill what 030 leaves uncovered |
| `validFootprints*.geojson` | `scanMasks` valid-data (mask 1..254) footprints of the partial frames; some are entirely mask 0 (partially focused) and cover nothing |
| `both/`, `descending/` | cycle 030 only (acquisition footprints) |
| `ascending/` | cycle 030 + 029/031 fill, valid footprints; built keeping the earlier ascending picks (`--keepFrom`), so only tiles with gaps changed |

`ascending/` was built with:

    globalGCOVTiles --catalogue cycle30/catalogue.geojson \
        --fillCatalogues cycle30/catalogue029.geojson cycle30/catalogue031.geojson \
        --validFootprints cycle30/validFootprints.geojson cycle30/validFootprints029.geojson cycle30/validFootprints031.geojson \
        --keepFrom <previous ascending> --direction A --allowV --exclude5MHzOcean --out cycle30/ascending
    coverageMaps cycle30/ascending

Rebuilt 2026-09-28 with the antimeridian fix (ASF writes footprints crossing 180 with lon > 180;
they were missing from the W180 tiles): 10 W180 tiles gained granules, 12 jobs.

Rebuilt 2026-09-28 (gap fill): the finished run's progress.png -> gapsFromPicture (land, ice
shelves, lakes, Caspian) -> ascending/fillGaps.geojson -> globalGCOVTiles --keepFrom <previous>
--fillGaps ... with valid footprints of ALL cycle-030/031 ascending granules (029 partly scanned;
ASF returned 503s). 211 tiles gained 581 granules (030: 107, 029: 405, 031: 69); 243 jobs, incl. the cap.

Uncovered land: 2.5% (cycle 030 only) -> 0.7% (with fill): the Antarctic pole hole beyond ascending
reach, a sliver near 85S 180, and 5 MHz stripes in the Sahara.
