# cycle30 tilings

| file | what |
|---|---|
| `catalogue.geojson` | cycle 030 GCOV granules (ASF search), `L2_PR`, one per scene |
| `catalogue029.geojson`, `catalogue031.geojson`, `catalogue027.geojson` | cycles 029, 031 and 027 (08-10; only 4247 granules so far), used only to fill what 030 leaves uncovered, in that order. All four catalogues are NISAR_L2_GCOV_PROVISIONAL_V1, CRID P05023 |
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

Rebuilt 2026-09-28 again (option B): (1) the northern Sahara stripe is ascending track 130, which is
not acquired north of ~22 N in cycles 029-031 but is in cycle 027 (frames 014-019, 08-10); gaps left
after the picture fill (fillGaps-left.geojson) filled with 027 added; (2) every land gap left by the
kept granules' valid data (all cycle-030/031 ascending granules scanned; fillGaps-valid.geojson),
e.g. datatake-break stripes in the central US and specks in Antarctica. 470 jobs change against the
tiling of 7c3546a; the cap does not change. Uncovered land 0.5%.

Uncovered land: 2.5% (cycle 030 only) -> 0.7% (with fill): the Antarctic pole hole beyond ascending
reach, a sliver near 85S 180, and 5 MHz stripes in the Sahara.
