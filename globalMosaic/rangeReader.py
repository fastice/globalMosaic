'''
Read chosen HDF5 datasets of a remote NISAR granule without downloading the rest of the file.

GDAL's /vsicurl reads ahead in large blocks, which is right for data stored contiguously but not
for a GCOV mask: its 4556 tiny chunks are scattered through the whole 7 GB file, so every mask
read pulled a 10 MB block of other layers and a 40 MHz granule cost 5.4 GB for 1.7 GB of HHHH.
Here h5py reads the file through RangeFile: before each strip the exact byte ranges of the chunks
it needs are fetched in parallel (prefetch), and h5py then decodes them from memory; metadata
comes in small blocks on demand.
'''
import io
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import requests

_local = threading.local()


def session():
    ''' One requests session per process (Earthdata login via ~/.netrc, cookies kept). '''
    if getattr(_local, 's', None) is None:
        s = requests.Session()
        a = requests.adapters.HTTPAdapter(pool_connections=16, pool_maxsize=16, max_retries=3)
        s.mount('https://', a)
        _local.s = s
    return _local.s


class RangeFile(io.RawIOBase):
    ''' Read-only, seekable file over HTTP range requests, for h5py. '''

    def __init__(self, url, block=256 * 1024, threads=8, piece=16 * 1024 * 1024, gap=64 * 1024):
        self.url, self.block, self.threads, self.piece, self.gap = url, block, threads, piece, gap
        self.pos = 0
        self.blocks = {}                    # metadata: block index -> bytes
        self.spans = []                     # prefetched (start, bytes), sorted by start
        self.downloaded = 0
        self._resolve()

    def _resolve(self):
        ''' Follow the Earthdata redirects once; range requests then go to the signed URL. '''
        r = session().get(self.url, headers={'Range': 'bytes=0-0'}, stream=True, timeout=120)
        r.raise_for_status()
        self.signed = r.url
        self.size = int(r.headers['Content-Range'].split('/')[-1])
        r.close()

    def _get(self, start, stop):
        for attempt in range(4):
            r = session().get(self.signed, headers={'Range': f'bytes={start}-{stop - 1}'}, timeout=300)
            if r.status_code in (401, 403) and attempt < 3:
                self._resolve()             # signed URL expired
                continue
            r.raise_for_status()
            self.downloaded += len(r.content)
            return r.content
        raise IOError(f'range {start}-{stop} of {self.url}')

    def prefetch(self, ranges):
        ''' Fetch these (offset, size) byte ranges (merged when close, split into pieces) in parallel. '''
        merged = []
        for o, n in sorted(ranges):
            if merged and o - merged[-1][1] <= self.gap and o + n - merged[-1][0] <= self.piece:
                merged[-1][1] = max(merged[-1][1], o + n)
            else:
                merged.append([o, o + n])
        with ThreadPoolExecutor(self.threads) as pool:
            data = list(pool.map(lambda ab: self._get(*ab), merged))
        self.spans = sorted(zip((a for a, _ in merged), data), key=lambda x: x[0])
        self._starts = [a for a, _ in self.spans]

    def drop(self):
        self.spans, self._starts = [], []

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        out = memoryview(b)
        got = 0
        while got < n:
            p = self.pos + got
            chunk = self._fromSpans(p, n - got)
            if chunk is None:
                i = p // self.block
                if i not in self.blocks:
                    self.blocks[i] = self._get(i * self.block, min((i + 1) * self.block, self.size))
                blk = self.blocks[i]
                chunk = blk[p - i * self.block:p - i * self.block + n - got]
            out[got:got + len(chunk)] = chunk
            got += len(chunk)
        self.pos += n
        return n

    def _fromSpans(self, p, want):
        import bisect
        j = bisect.bisect_right(getattr(self, '_starts', []), p) - 1
        if j < 0:
            return None
        a, d = self.spans[j]
        if p >= a + len(d):
            return None
        return d[p - a:p - a + want]


def chunkIndex(ds):
    ''' {(row0, col0): (byte offset, size)} of a chunked dataset's stored chunks. '''
    out = {}
    for i in range(ds.id.get_num_chunks()):
        c = ds.id.get_chunk_info(i)
        out[tuple(c.chunk_offset)] = (c.byte_offset, c.size)
    return out


def rangesFor(index, chunks, r0, r1, c0, c1):
    ''' Byte ranges of the chunks that overlap rows r0:r1, columns c0:c1. '''
    ch, cw = chunks
    out = []
    for rr in range(r0 // ch * ch, r1, ch):
        for cc in range(c0 // cw * cw, c1, cw):
            v = index.get((rr, cc))
            if v is not None:
                out.append(v)
    return out
