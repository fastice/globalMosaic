#!/usr/bin/env bash
# Build fp16read against whatever GDAL+HDF5 stacks are present. -Wl,-rpath is what lets a
# conda-linked binary run without LD_LIBRARY_PATH, which otherwise silently shadows the
# system libraries for everything else in the shell.
set -eu
cd "$(dirname "$0")"
built=0
if [ -f /usr/include/gdal/gdal.h ] && [ -f /usr/include/hdf5/serial/hdf5.h ]; then
    gcc -O3 -o fp16read.sys fp16read.c -I/usr/include/gdal -I/usr/include/hdf5/serial \
        -L/usr/lib/$(uname -m)-linux-gnu -lgdal -lhdf5_serial -lm && \
        { echo "built fp16read.sys (system stack)"; built=1; }
fi
P="${CONDA_PREFIX:-/opt/conda}"
if [ -f "$P/include/gdal.h" ] && [ -f "$P/include/hdf5.h" ]; then
    gcc -O3 -o fp16read.conda fp16read.c -I"$P/include" \
        -L"$P/lib" -Wl,-rpath,"$P/lib" -lgdal -lhdf5 -lm && \
        { echo "built fp16read.conda (conda stack: $P)"; built=1; }
fi
[ "$built" = 1 ] || { echo "no usable GDAL+HDF5 stack found" >&2; exit 1; }
