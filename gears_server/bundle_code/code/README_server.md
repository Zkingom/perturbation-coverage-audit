# GEARS robustness run on a rented GPU server (RPE1 essential, 21 frozen splits)

Protocol: v0.6 section 10 (inherited by v0.7). GEARS 0.1.2, hidden_size 64, 20 epochs, uncertainty off, best epoch by
validation top-20 DE MSE; training data capped at 16 cells per perturbation (fixed seed), all control cells.

## 0. Upload and unpack (on the server)
    mkdir -p ~/gears && cd ~/gears
    tar -xzf gears_bundle_v07.tar.gz          # creates ./gears_bundle_v07
    cd gears_bundle_v07 && sha256sum -c SHA256SUMS --quiet && echo BUNDLE OK

## 1. Environment + self-check (once; ~10-20 min)
    bash code/setup_gears_env.sh               # optional argument: a pip mirror URL
    # must end with "SETUP PASS"; details in logs/setup_selfcheck.json

## 2. Preprocessing (once; GEARS DE genes + per-cell graphs)
    source venv_gears/bin/activate
    nohup python code/gears_server_run.py prepare > logs/prepare.txt 2>&1 &
    tail -f logs/prepare.txt                   # ends with "prepare done {...}"

## 3. Smoke test: one split
    nohup python code/gears_server_run.py run --splits 0 > logs/smoke.txt 2>&1 &
    # ends with "split 00 done {...}": send me that line (time, peak RSS, peak GPU memory)

## 4. All 21 splits (only after the smoke test was reviewed)
    nohup python code/gears_server_run.py run > logs/run_all.txt 2>&1 &
    # interrupted? run the same command again: verified splits are skipped, an unfinished split restarts

## 5. Download
    python code/gears_server_run.py verify     # writes results/SHA256SUMS
    tar -czf gears_results_v07.tar.gz results logs
    # download gears_results_v07.tar.gz; it is re-hashed locally before use

Use `tmux` or `nohup` so the run survives a dropped SSH connection. Do not edit files in data/, priors/, splits/.
