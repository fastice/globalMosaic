# cycle30gunw: GUNW (interferogram) tilings for coherence mosaics

| file | what |
|---|---|
| `catalogue.geojson` | reference cycle 030, 12-day pairs (030 -> 031) |
| `catalogue029/028/027.geojson` | 12-day pairs with reference cycle 029 / 028 / 027: fill, in that order |
| `catalogue*.all.geojson` | everything ASF lists for that reference cycle (incl. 24/36-day pairs) |
| `ascending/`, `descending/` | tilings, direction by file name (not ASF's flightDirection) |

Built 2026-09-29 with `globalGCOVTiles --product GUNW` searches, filtered to secondary cycle =
reference + 1, then `globalGCOVTiles --catalogue catalogue.geojson --fillCatalogues catalogue029
catalogue028 catalogue027 --direction A|D`. Cycle 031 is still being acquired, so 030 -> 031 pairs
are incomplete and the fill cycles carry more. Mosaic: `mosaicCoherence` / `test/runCoherence.sh`
(unwrappedInterferogram/HH/coherenceMagnitude, 80 m, averaged where pairs overlap).
