# NISAR L-band backscatter, global ascending mosaic (Google Earth)

Radiometrically terrain-corrected backscatter (gamma0) from NISAR L2 GCOV products, mosaicked
globally from **ascending passes only** and packaged as a Google Earth tiled overlay.

## The data

| | |
|---|---|
| Sensor / product | NISAR L-band SAR, L2 GCOV (provisional, CRID P05023) |
| Quantity | gamma0 (RTC), shown in dB |
| Polarization | HH; VV where a granule has no H channel (e.g. some 5 MHz modes) |
| Passes | ascending only |
| Time | cycle 030 (September 2026), with gaps filled from cycles 023-031 (mid-June to late September 2026) |
| Mosaic spacing | 2.4 arcsec lat/lon (~75 m); 75 m polar stereographic south of 84 S |
| Display stretch | grey scale, -24 dB (black) to -1 dB (white); no data is transparent |

Land that has no ascending acquisition in any of these cycles (e.g. parts of northern
Thailand/Laos, stripes in the Sahara, the western Aleutians, Faroe, Shetland, and the Antarctic
pole hole south of 87 S) is left empty; see `notAcquiredAscending/` for the full list.

## Files

| File | Contents | Google Earth zoom | Finest pixel |
|---|---|---|---|
| `ge_base.tar` | `googleEarth/doc.kml` (top level) + `googleEarth/base/` | 0-8 | ~300 m |
| `ge_land.tar` | `googleEarth/land/` (land only, lightly smoothed) | 4-9 (images at 9) | ~150 m |
| `ge_detail.tar` | `googleEarth/detail/` (40/77 MHz land, ice-sheet margins, glaciers) | 6-10 (images at 10) | ~75 m |
| `ge_detail20.tar` (optional) | `googleEarth/detail20/` (the rest of the land covered by 20 MHz data; not ocean or sea ice, not 5 MHz-only land, not the ice-sheet interiors) | 6-10 (images at 10) | ~75 m |

Tiles are 512 x 512 PNGs in a KML superoverlay. Google Earth loads them as you zoom.

## Opening it

1. Unpack the tars **into the same folder** (each carries the same top `doc.kml`):

       tar -xf ge_base.tar
       tar -xf ge_land.tar
       tar -xf ge_detail.tar
       tar -xf ge_detail20.tar      # optional

2. In Google Earth Pro: **File > Open** and pick one of:
   - `googleEarth/doc.kml`: the full product (all layers)
   - `googleEarth/base/doc.kml`: the low-resolution global layer only

Base works on its own and with any combination of the others. Land, detail and detail20 hold images
only at their finest level (their coarser levels are just links), so on their own they show nothing
until you zoom in that far; base supplies the coarser zooms.
If one of them is missing, the top `doc.kml` shows that entry as a broken link but displays the
rest normally.

If the tars were downloaded in pieces (`*.tar.part-00`, `-01`, ...), join them first:

    cat ge_detail.tar.part-* > ge_detail.tar
    md5sum -c ge_detail.tar.md5

## Notes

- Over the ocean (sea ice) the tiles sit 10 m above sea level so that Google Earth's water surface
  does not hide them; coastal tiles are split into a land part (draped on the terrain) and a sea
  part. Leave **Terrain** on; it looks correct either way.
- Built with the globalMosaic package (github.com/fastice/globalMosaic): `runGeomosaicTiles`
  (mosaic) and `makeGoogleEarthTiered` (Google Earth layers); tiling and fill history in
  `cycle30/README.md`.
