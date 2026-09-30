# Land not acquired ascending (cycles 023-031)

For the acquisition planners: land (Natural Earth 10 m + Antarctic ice shelves, south of 77.5 N)
that no ascending GCOV granule of cycles 023-031 (2026-06-17 to cycle 031) covers, by acquisition
footprint. These gaps cannot be filled from the archive.

    notAcquired --catalogues cycle30/catalogue0{23..29}.geojson cycle30/catalogue.geojson cycle30/catalogue031.geojson \
        --direction A --out cycle30/notAcquiredAscending

| file | what |
|---|---|
| `notAcquired.csv` | one row per connected area >= 4 km2 (2 km cells), largest first: area, country, nearest named place, representative point, bounds, share covered by descending passes |
| `notAcquired.geojson` | the same areas as polygons |
| `notAcquired.png` | map |

298 areas, 442,466 km2. 251,000 km2 is the Antarctic pole hole beyond ascending reach (expected).
The rest: northern Thailand/Laos (37,000 km2), coastal Ghana, many stripes between tracks in the
Sahara (Algeria, Niger, Mali, Mauritania, Chad, Libya, Sudan, Egypt), Reunion and Mauritius (no
descending either), the western Aleutians, Faroe, Shetland, Aldabra, the Maldives and many
Pacific islands.

Acquired but with no valid data there (partially focused, mask 0, in every candidate): Costa Brava
(Spain), Ishigaki, Miyako and Okinawa (Japan). These are not in the list above, which uses
acquisition footprints.
