#!/usr/bin/env bash
# run_all.sh v3 — 冻结包 code/ 一个字节都不改；内存/速度补丁只在运行时通过 PYTHONPATH 下的 sitecustomize 注入。
#   每个 split：轻量进程预算网络(precompute_nets.py) -> 正式运行 code/gears_server_run.py run --splits i
#   参数(batch size / epochs / cell cap / seed / hidden / uncertainty)全部来自冻结代码本身，这里不传任何参数；串行运行。
#   出错不自动修复：失败的 split 记录在 logs/run_all_timing.txt，继续下一个；全部结束后按冻结方案 verify + 打包。
set -u
B=/root/autodl-tmp/gears/gears_bundle_v07
RT=/root/autodl-tmp/gears/runtime_patches
PY=$B/venv_gears/bin/python
START=${START:-0}
cd "$B" || { echo "找不到 $B"; exit 1; }
sha256sum -c SHA256SUMS --quiet || { echo "BUNDLE CHECKSUM FAILED" | tee /root/autodl-tmp/gears/ERROR.txt; exit 1; }
mkdir -p results logs
one() {
  i=$1; s=$(printf "%02d" "$i")
  if [ -f results/split_$s/done.json ]; then echo "$(date '+%F %T') split $s 已完成，跳过"; return 0; fi
  PYTHONPATH=$RT "$PY" $RT/precompute_nets.py "$i" > logs/pre_$s.txt 2>&1 || echo "$(date '+%F %T') split $s precompute failed" | tee -a logs/run_all_timing.txt
  echo "$(date '+%F %T') split $s 开始"
  t0=$(date +%s)
  PYTHONPATH=$RT "$PY" code/gears_server_run.py run --splits "$i" > logs/split_$s.txt 2>&1
  rc=$?
  echo "$(date '+%F %T') split $s 结束 rc=$rc 用时 $(( ($(date +%s)-t0)/60 )) 分钟" | tee -a logs/run_all_timing.txt
  ( flock 9; tar czf /root/autodl-tmp/gears/RESULTS_latest.tmp.tar.gz results logs && mv /root/autodl-tmp/gears/RESULTS_latest.tmp.tar.gz /root/autodl-tmp/gears/RESULTS_latest.tar.gz ) 9>/root/autodl-tmp/gears/.pack.lock
}
for i in $(seq $START 20); do one $i; done
echo "$(date '+%F %T') 全部 split 处理完毕"
sha256sum -c SHA256SUMS --quiet && echo "BUNDLE STILL OK"
"$PY" code/gears_server_run.py verify > logs/verify.txt 2>&1; echo "verify rc=$?"
tar -czf gears_results_v07.tar.gz results logs && sha256sum gears_results_v07.tar.gz > gears_results_v07.tar.gz.sha256 && echo "PACKED gears_results_v07.tar.gz"
