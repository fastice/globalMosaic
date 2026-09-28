#!/usr/bin/env python3
"""
Per-interface network rate, like nload, with no dependencies and no root.

Reads /proc/net/dev, which every container exposes. Use it during a tile run to tell whether
you are network-bound or CPU-bound: if the rate plateaus while cores sit idle, add processes;
if cores are pinned and the rate is below the instance's ceiling, you are CPU-bound.

    python netmon.py                 # all busy interfaces, 2 s samples
    python netmon.py -i eth0 -n 30   # one interface, 30 samples then stop
"""
import argparse
import sys
import time


def readDev():
    stats = {}
    with open('/proc/net/dev') as f:
        for line in f.readlines()[2:]:
            name, _, rest = line.partition(':')
            v = rest.split()
            stats[name.strip()] = (int(v[0]), int(v[8]))   # rx bytes, tx bytes
    return stats


def human(bps):
    for unit in ('B/s', 'KB/s', 'MB/s', 'GB/s'):
        if bps < 1024 or unit == 'GB/s':
            return '%7.2f %s' % (bps, unit)
        bps /= 1024.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--interface', default=None, help='one interface [all busy ones]')
    ap.add_argument('-t', '--interval', type=float, default=2.0, help='sample seconds [2]')
    ap.add_argument('-n', '--count', type=int, default=0, help='stop after N samples [forever]')
    args = ap.parse_args()

    prev = readDev()
    t0 = time.time()
    peak = {}
    total = {}
    k = 0
    try:
        while args.count == 0 or k < args.count:
            time.sleep(args.interval)
            now = readDev()
            dt = time.time() - t0
            t0 = time.time()
            out = []
            for name in sorted(now):
                if args.interface and name != args.interface:
                    continue
                if name == 'lo' and not args.interface:
                    continue
                rx = (now[name][0] - prev[name][0]) / dt
                tx = (now[name][1] - prev[name][1]) / dt
                if not args.interface and rx < 1024 and tx < 1024:
                    continue                      # skip idle interfaces
                peak[name] = max(peak.get(name, 0), rx)
                total[name] = total.get(name, 0) + (now[name][0] - prev[name][0])
                out.append('%s rx %s  tx %s  (peak rx %s)'
                           % (name, human(rx), human(tx), human(peak[name])))
            print(' | '.join(out) if out else '(idle)', flush=True)
            prev = now
            k += 1
    except KeyboardInterrupt:
        pass
    for name in sorted(total):
        print('%s: %.2f GB received' % (name, total[name] / 1024.0 ** 3), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
