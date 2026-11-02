# NISAR dual-pol false-colour composite

A three-band RGB composite of like-pol and cross-pol backscatter and their ratio, following
Earth Big Data's openSEPPO `-dpratio` recipe (`doc/ratio.md` in `EarthBigData/openSEPPO`).
Built by `ratioMosaic --mode rgb`, usually through `test/runRatioTiles.sh`.

## What the bands are

| band | content |
|---|---|
| R | like-pol gamma0, as DN |
| G | cross-pol gamma0, as DN |
| B | the like/cross ratio, `DN(like) - DN(cross) + 127` |

DN is a linear dB scale, openSEPPO's `power_to_dn_uint8`:

    dB = -31.15 + DN * 0.15          -31.0 dB -> DN 1, +7.1 dB -> DN 255, DN 0 = no data

Byte, **no data 0**, with a per-band display stretch applied when the tiles were built (see
below). The ratio is computed from the *clamped* DNs, as openSEPPO does, not by scaling a dB
difference -- cross-pol below -31 dB pins to DN 1, which is routine over desert and water, and
the two disagree exactly there.

## Reading it

- **Pale yellow** -- vegetation. Strong volume scattering from canopy, so cross-pol is high and
  the ratio low.
- **Purple / indigo** -- bare or arid ground. Little depolarisation, cross-pol near the noise
  floor, ratio high.
- **Grey** -- single-pol acquisitions, where there is no cross-pol and so no ratio. All three
  bands carry the like-pol backscatter instead, so the ground still shows. This is how the ice
  sheets appear, and the ~16% of the globe acquired as 5 MHz VV.
- **Green-yellow patches**, notably over western Russia and Kazakhstan -- **RFI**, not land
  cover. RFI inflates cross-pol, which is the only way G exceeds B.

Blue never falls below mid scale before stretching: cross-pol is always weaker than like-pol,
so B is confined to [127, 255] by construction. Anything below 127 means cross-pol exceeded
like-pol, which surface scattering cannot do.

## What is excluded

- **The ratio is dropped where cross-pol is at or below -28 dB** (`--minCross`). NISAR's noise
  floor is about -25 dB; below it the ratio is a lower bound on a quantity that was not
  measured, and open water would otherwise be the brightest thing in the image. Those pixels
  fall back to grey, so backscatter is still shown where the ratio is not.
- **Single-pol granules are excluded from the cross-pol cover entirely** (`easeTiles`/
  `globalGCOVTiles --crossPol`), because they have no cross-pol term. About 70% of NISAR
  granules are dual-pol, so the composite covers correspondingly less ground than a plain
  backscatter mosaic.
- **Polar caps are a separate CRS.** The lat/lon tiling puts everything beyond 84 degrees in
  polar stereographic tiles (EPSG 3031 south, 3413 north), which are not in the lat/lon VRT.

## The stretch

There is no published scaling to match: the Earth Big Data map credits "openSEPPO v0.9.0 and
SEPPO v7.13.0 for scaling", and SEPPO's scaling is not documented. The numbers below were
settled by eye and each does something specific; they are a look, not a standard.

Limits are derived ONCE and applied to every tile with `--stretchLimits`, never per tile --
per-tile percentiles give each tile its own and the mosaic comes out patchwork. Colour limits
come from dual-pol LAND only (60S-70N, grey pixels excluded), grey limits from the grey
population separately, sampling with `-r nearest` rather than `average` so the limits match the
full-resolution speckle the tiles actually contain.

For the cycle-30 ascending 2.4 arcsec product:

    --minCross -28 --stretchLimits 65 185 28 155 122 214 --greyLimits 20 205

## Same-granule caveat

A ratio must divide one acquisition by itself. That holds when both passes draw from one cover,
which is what `--ratioCover` (like-pol) and `--crossPol` (cross-pol) are for. If the like-pol
mosaic was built on an unrestricted cover instead -- as the cycle-30 z11 product was -- then
wherever that cover picked a single-pol granule the cross-pol pass picked a different one, and
those pixels compare two acquisitions. In practice single-pol is almost entirely ice sheets,
which have no dual-pol alternative to disagree with, so the affected fraction is small; but it
is not zero and it is not flagged in the data.
