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

Uncovered land: 2.5% (cycle 030 only) -> 0.7% (with fill): the Antarctic pole hole beyond ascending
reach, a sliver near 85S 180, and 5 MHz stripes in the Sahara.
