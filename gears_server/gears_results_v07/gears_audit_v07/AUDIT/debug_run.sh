#!/usr/bin/env bash
pkill -f run_all.sh; pkill -f gears_server_run.py; pkill -f precompute_nets.py; sleep 3
B=/root/autodl-tmp/gears/gears_bundle_v07; D=/root/autodl-tmp/gears/debug_bundle
rm -rf $D; mkdir -p $D/code $D/results $D/logs
cp $B/code/gears_server_run.py $D/code/
for x in data priors splits venv_gears work; do ln -s $B/$x $D/$x; done
cd $D
( while true; do echo "$(date +%T) $(( $(grep -E '^total_rss ' /sys/fs/cgroup/memory/memory.stat | cut -d' ' -f2)/1048576 ))MB" >> /root/autodl-tmp/gears/debug_mem.txt; sleep 1; done ) &
S=$!
DEBUG_EPOCHS=1 PYTHONPATH=/root/autodl-tmp/gears/debug_patches $B/venv_gears/bin/python code/gears_server_run.py run --splits 0 > /root/autodl-tmp/gears/debug_run.txt 2>&1
echo "EXIT rc=$?" >> /root/autodl-tmp/gears/debug_run.txt
kill $S
