/*
  fp16read -- time the three ways C can get a float16 or float32 raster out of HDF5.

  This mirrors what geomosaic actually does (geoMosaic/gcovMosaic.c): GDAL reports a float16
  band as Float32 and asks libhdf5 to convert, which lands in a generic software path; reading
  the raw 16-bit values and expanding them through a 65536-entry lookup table is much faster.
  The table is used instead of F16C intrinsics so no architecture-specific compiler flag is
  needed and the code still builds on arm64.

    fp16read gdal  <file> <dataset> <x> <y> <w> <h>   read via GDAL RasterIO
    fp16read h5    <file> <dataset> <x> <y> <w> <h>   read via HDF5, converting in HDF5
    fp16read lut   <file> <dataset> <x> <y> <w> <h>   read raw fp16 bits, expand with the LUT
    fp16read write <file> <dataset> <w> <h> <16|32> <deflate>   write a test dataset

  Prints: seconds, then Mpx/s, on one line.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <hdf5.h>
#include "gdal.h"
#include "cpl_conv.h"

static float f16Table[65536];

static void buildF16Table(void)
{
	unsigned int h, sign, exp, mant, bits;
	for (h = 0; h < 65536; h++)
	{
		sign = (h & 0x8000u) << 16;
		exp = (h >> 10) & 0x1fu;
		mant = h & 0x3ffu;
		if (exp == 0)
		{
			if (mant == 0)
			{
				bits = sign;
			}
			else
			{
				exp = 127 - 15 + 1;
				while ((mant & 0x400u) == 0)
				{
					mant <<= 1;
					exp--;
				}
				mant &= 0x3ffu;
				bits = sign | (exp << 23) | (mant << 13);
			}
		}
		else if (exp == 31)
		{
			bits = sign | 0x7f800000u | (mant << 13);
		}
		else
		{
			bits = sign | ((exp - 15 + 127) << 23) | (mant << 13);
		}
		memcpy(&f16Table[h], &bits, 4);
	}
}

static double now(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

static void report(double dt, long npx)
{
	printf("%.4f %.2f\n", dt, npx / 1e6 / dt);
}

int main(int argc, char **argv)
{
	char path[4096];
	double t0;
	long npx;

	if (argc < 2)
	{
		fprintf(stderr, "usage: see the comment at the top of fp16read.c\n");
		return 2;
	}
	if (strcmp(argv[1], "write") == 0)
	{
		/* write <file> <dataset> <w> <h> <bits> <deflate> */
		hid_t f, s, d, p, mt, ft;
		hsize_t dims[2], ch[2];
		int w = atoi(argv[4]), h = atoi(argv[5]), bits = atoi(argv[6]), dfl = atoi(argv[7]);
		long i;
		float *buf = (float *)malloc(sizeof(float) * (size_t)w);
		dims[0] = h;
		dims[1] = w;
		ch[0] = (hsize_t)(h < 512 ? h : 512);
		ch[1] = (hsize_t)(w < 512 ? w : 512);
		/* realistic-ish backscatter: log-normal, so deflate behaves like it does on real data */
		for (i = 0; i < w; i++)
		{
			buf[i] = (float)exp(-3.0 + 1.5 * ((double)rand() / RAND_MAX - 0.5) * 4.0);
		}
		/* H5T_IEEE_F16LE only exists from HDF5 1.14.4; both stacks here are older (1.10.10
		   and 1.14.3), and the slim GCOV products predate it too. Build IEEE binary16 by hand
		   instead, which works on every version and is what lets 1.10 read those products. */
		if (bits == 16)
		{
			ft = H5Tcopy(H5T_IEEE_F32LE);
			H5Tset_size(ft, 2);
			H5Tset_precision(ft, 16);
			H5Tset_fields(ft, 15, 10, 5, 0, 10);
			H5Tset_ebias(ft, 15);
		}
		else
		{
			ft = H5Tcopy(H5T_IEEE_F32LE);
		}
		mt = H5T_NATIVE_FLOAT;
		f = H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		s = H5Screate_simple(2, dims, NULL);
		p = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(p, 2, ch);
		if (dfl > 0)
		{
			H5Pset_deflate(p, dfl);
		}
		d = H5Dcreate2(f, argv[3], ft, s, H5P_DEFAULT, p, H5P_DEFAULT);
		t0 = now();
		{
			hsize_t off[2], cnt[2];
			hid_t ms, fs;
			cnt[0] = 1;
			cnt[1] = (hsize_t)w;
			ms = H5Screate_simple(2, cnt, NULL);
			fs = H5Dget_space(d);
			for (i = 0; i < h; i++)
			{
				off[0] = (hsize_t)i;
				off[1] = 0;
				H5Sselect_hyperslab(fs, H5S_SELECT_SET, off, NULL, cnt, NULL);
				H5Dwrite(d, mt, ms, fs, H5P_DEFAULT, buf);
			}
			H5Sclose(ms);
			H5Sclose(fs);
		}
		H5Fflush(f, H5F_SCOPE_GLOBAL);
		npx = (long)w * h;
		report(now() - t0, npx);
		H5Tclose(ft);
		H5Dclose(d);
		H5Pclose(p);
		H5Sclose(s);
		H5Fclose(f);
		free(buf);
		return 0;
	}

	/* read modes: <file> <dataset> <x> <y> <w> <h> */
	{
		int x = atoi(argv[4]), y = atoi(argv[5]), w = atoi(argv[6]), h = atoi(argv[7]);
		npx = (long)w * h;
		if (strcmp(argv[1], "gdal") == 0)
		{
			GDALDatasetH ds;
			GDALRasterBandH b;
			float *buf = (float *)malloc(sizeof(float) * (size_t)npx);
			GDALAllRegister();
			snprintf(path, sizeof(path), "HDF5:\"%s\"://%s", argv[2], argv[3]);
			ds = GDALOpen(path, GA_ReadOnly);
			if (ds == NULL)
			{
				fprintf(stderr, "cannot open %s\n", path);
				return 1;
			}
			b = GDALGetRasterBand(ds, 1);
			t0 = now();
			GDALRasterIO(b, GF_Read, x, y, w, h, buf, w, h, GDT_Float32, 0, 0);
			report(now() - t0, npx);
			GDALClose(ds);
			free(buf);
			return 0;
		}
		else
		{
			hid_t f, d, fs, ms, ft;
			hsize_t off[2], cnt[2];
			int lut = (strcmp(argv[1], "lut") == 0);
			float *out = (float *)malloc(sizeof(float) * (size_t)npx);
			unsigned short *raw = NULL;
			long i;
			f = H5Fopen(argv[2], H5F_ACC_RDONLY, H5P_DEFAULT);
			if (f < 0)
			{
				fprintf(stderr, "cannot H5Fopen %s\n", argv[2]);
				return 1;
			}
			d = H5Dopen2(f, argv[3], H5P_DEFAULT);
			ft = H5Dget_type(d);
			off[0] = (hsize_t)y;
			off[1] = (hsize_t)x;
			cnt[0] = (hsize_t)h;
			cnt[1] = (hsize_t)w;
			fs = H5Dget_space(d);
			H5Sselect_hyperslab(fs, H5S_SELECT_SET, off, NULL, cnt, NULL);
			ms = H5Screate_simple(2, cnt, NULL);
			if (lut && H5Tget_size(ft) == 2)
			{
				buildF16Table();
				raw = (unsigned short *)malloc(sizeof(unsigned short) * (size_t)npx);
				t0 = now();
				/* read using the FILE's own type: HDF5 copies bits, it does not convert */
				H5Dread(d, ft, ms, fs, H5P_DEFAULT, raw);
				for (i = 0; i < npx; i++)
				{
					out[i] = f16Table[raw[i]];
				}
				report(now() - t0, npx);
				free(raw);
			}
			else
			{
				t0 = now();
				H5Dread(d, H5T_NATIVE_FLOAT, ms, fs, H5P_DEFAULT, out);
				report(now() - t0, npx);
			}
			H5Sclose(ms);
			H5Sclose(fs);
			H5Tclose(ft);
			H5Dclose(d);
			H5Fclose(f);
			free(out);
			return 0;
		}
	}
}
