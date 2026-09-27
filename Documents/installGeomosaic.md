# installGeomosaic

Clones (or pulls) the GrIMP C repositories that `geomosaic` needs, builds them, installs the
binary, and then **verifies the binary can actually do what globalMosaic requires**.

`geomosaic` is the only GrIMP binary globalMosaic drives. GIT64 is not one repository — each
component is its own — so five must be checked out side by side under a single root, because
`mosaicSource/Makefile` reaches across to `$(PROGDIR)/<repo>` (see `GEOMOSAICDIRS`):

| repo | supplies |
|---|---|
| `mosaicSource` | `geoMosaic/`, `common/`, `landsatMosaic/` — the binary itself |
| `gdalIO` | GDAL / GeoTIFF I/O |
| `clib` | byte-swapped flat-file I/O |
| `cRecipes` | SVD, interpolation |
| `triangle` | Shewchuk triangulation |

All five are public on `github.com/fastice`, so a fresh instance needs no SSH key.

## Usage

```bash
installGeomosaic                       # console script, after pip install -e .
install/installGeomosaic               # or run it directly, no Python needed
```

| option | meaning | default |
|---|---|---|
| `--gitroot DIR` | where the repos live | `$HOME/progs/GIT64` |
| `--bindir DIR` | where `geomosaic` is installed | `$HOME/bin/$(uname -m)` |
| `--branch NAME` | force one branch for every repo | per-repo table |
| `--ssh` | clone over ssh instead of https | https |
| `--jobs N` | parallel make jobs | `nproc` |
| `--no-pull` | build what is on disk, do not fetch | pull |
| `--force` | discard local modifications when pulling | never discard |
| `--check-only` | report what would happen, build nothing | — |

`GEOMOSAIC_GITROOT` and `GEOMOSAIC_BINDIR` work as environment equivalents.

Typical fresh instance:

```bash
sudo apt-get install -y build-essential git libgdal-dev libhdf5-dev libproj-dev
installGeomosaic
export PATH="$HOME/bin/$(uname -m):$PATH"
export OPENBLAS_NUM_THREADS=1
```

## Branch pinning — read this before it bites you

`clib`, `cRecipes` and `triangle` track `master`. **`gdalIO` and `mosaicSource` are pinned to
`epsg-proj-support`**, because that is where the arbitrary-EPSG work lives — including
`-epsg 4326` (geographic output) and `dem none`, both of which globalMosaic cannot run without.

Two traps live here:

- The repos' default branches **differ**: `mosaicSource` uses `master`, `gdalIO` uses `main`.
  Do not assume one name.
- As of 2026-09-27 **`epsg-proj-support` is not published** on either remote, so a clone from
  GitHub cannot build a usable `geomosaic` yet. The script detects this with `git ls-remote`
  and stops with a message naming the repo, rather than silently falling back to the default
  branch and building a `geomosaic` that lacks the features.

Once that branch is merged, set both entries in the `BRANCH` table to each repo's default.

## The capability check is the real guarantee

A branch name can be stale, a merge can drop a feature, an old binary can sit earlier on
`PATH` — a successful `make` proves none of that. After installing, the script runs the binary
it just produced against a throwaway input with `-epsg 4326` and `dem none`, and requires the
line `grid units are DEGREES` in the output.

`geomosaic` is **expected to exit non-zero** on that check — the yaml lists no granules, so it
stops with `No range/Doppler or GCOV inputs`. What is being tested is that it got far enough
to accept `dem none` and resolve EPSG:4326; an older binary dies earlier trying to open
`none.geodat`. (The check writes to a file rather than a pipe, because `set -o pipefail` would
otherwise read that expected non-zero exit as a failed check.)

A failed capability check exits non-zero, so `subprocess.run(..., check=True)` catches it.

## Notes

- **Never discards your work.** A repo with modified tracked files is left alone with a
  warning, and pulls are `--ff-only`. `--force` overrides both.
- **System GDAL, not conda.** `geomosaic` builds against the system `libgdal`. If a conda lib
  directory is on `LD_LIBRARY_PATH` the script warns, since a different `libgdal` may be
  loaded at run time than it was built against.
- **`OPENBLAS_NUM_THREADS=1` matters at run time.** Every binary linking `libgdal` pulls in
  `libopenblas`, whose constructor spawns an idle thread pool sized to the machine before
  `main()` runs. No GrIMP code uses BLAS, and it can only be capped from the environment.
  Never use `OMP_NUM_THREADS=1` for this — it would serialise the genuinely parallel binaries.
