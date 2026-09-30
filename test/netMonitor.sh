#!/bin/bash
# Network receive rate, every few seconds, until Ctrl-C.
#   bash test/netMonitor.sh            # busiest interface, every 5 s
#   bash test/netMonitor.sh 10 eth0    # every 10 s, this interface
SEC=${1:-5}
IF=${2:-$(awk 'NR > 2 && $1 != "lo:" {gsub(":", "", $1); print $2, $1}' /proc/net/dev | sort -n | tail -1 | cut -d' ' -f2)}
rx() { awk -v i="$IF:" '$1 == i {print $2}' /proc/net/dev; }
echo "interface $IF, every $SEC s (Ctrl-C to stop)"
a=$(rx)
while sleep $SEC; do
    b=$(rx)
    echo "$(date +%T)  $(( (b - a) / SEC / 1000000 )) MB/s in"
    a=$b
done
