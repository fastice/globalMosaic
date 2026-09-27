#!/bin/bash
# Check that this machine can read NISAR GCOV granules, without geomosaic: Earthdata login,
# ASF S3 keys, a direct S3 read (AWS us-west-2 only) and an https (vsicurl) read.
#   bash test/checkAccess.sh
cd "$(dirname "$0")/.."
export PATH=$HOME/bin/$(uname -m):$PATH
echo "geomosaic: $(command -v geomosaic || echo NOT FOUND)"
python3 - <<'PY'
import netrc, os
from osgeo import gdal
from globalMosaic.runGeomosaicTiles import S3Keys, granulePath

print('GDAL', gdal.__version__)
drv = gdal.GetDriverByName('HDF5')
print('GDAL HDF5 driver:', 'yes' if drv else
      'NO -- geomosaic cannot read GCOV; on conda: conda install -c conda-forge libgdal-hdf5')
try:
    ok = netrc.netrc().authenticators('urs.earthdata.nasa.gov') is not None
    print('~/.netrc Earthdata entry:', 'yes' if ok else 'NO')
except Exception as e:
    print('~/.netrc:', e)

name = 'NISAR_L2_PR_GCOV_030_052_D_067_2005_DHDH_A_20260909T193347_20260909T193422_P05023_N_F_J_001'
url = ('https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L2_GCOV_PROVISIONAL_V1/'
       f'{name}/{name}.h5')


def tryRead(label, path):
    gdal.ErrorReset()
    f = gdal.VSIFOpenL(path, 'rb')
    if f is None:
        print(f'{label}: FAILED -- {gdal.GetLastErrorMsg()[:400]}')
        return
    data = gdal.VSIFReadL(1, 8, f)
    gdal.VSIFCloseL(f)
    print(f'{label}: OK (read {len(data)} bytes, HDF5 signature {data[1:4] == b"HDF"})')
    # what geomosaic does: GDALOpen of the file, then an HDF5 subdataset
    gdal.ErrorReset()
    try:
        ds = gdal.Open(f'HDF5:"{path}"://science/LSAR/GCOV/grids/frequencyA/HHHH')
        print(f'{label} GDAL HDF5 open: OK ({ds.RasterXSize} x {ds.RasterYSize})')
    except Exception as e:
        print(f'{label} GDAL HDF5 open: FAILED -- {str(e)[:300]}')


gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'EMPTY_DIR')
try:
    keys = S3Keys().env()
    print('ASF S3 keys: OK')
    for k, v in keys.items():
        gdal.SetConfigOption(k, v)
    tryRead('S3 read (granules: s3)', granulePath('s3', name, url))
    for k in keys:
        gdal.SetConfigOption(k, None)
except Exception as e:
    print(f'ASF S3 keys: FAILED -- {e}')

gdal.SetConfigOption('GDAL_HTTP_NETRC', 'YES')
gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE', '/tmp/checkAccess.cookies')
gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', '/tmp/checkAccess.cookies')
tryRead('https read (granules: vsicurl)', granulePath('vsicurl', name, url))
PY
