#!/bin/bash
# One-time setup: install globalMosaic, build geomosaic, check Earthdata credentials.
set -euo pipefail
cd "$(dirname "$0")/.."
pip install -e .
if ! python3 -c "from osgeo import gdal; import sys; sys.exit(gdal.GetDriverByName('HDF5') is None)"; then
    echo "GDAL has no HDF5 driver (geomosaic cannot read GCOV); installing the conda-forge plugin"
    conda install -y -c conda-forge libgdal-hdf5
fi
installGeomosaic
BIN=$HOME/bin/$(uname -m)
grep -qs "$BIN" ~/.bashrc || echo "export PATH=$BIN:\$PATH" >> ~/.bashrc
if ! grep -qs urs.earthdata.nasa.gov ~/.netrc; then
    echo "Add Earthdata credentials, then chmod 600 ~/.netrc:"
    echo "  machine urs.earthdata.nasa.gov login USER password PASS"
    exit 1
fi
echo "setup done; open a new shell (or: export PATH=$BIN:\$PATH)"
