# run_gwps_c0.ps1 -- user-run launcher for the GWPS C0 (protocol v0.7 section 3). Runs outside the sandbox because
# the single-cell file is on H: (exFAT, not reachable from the sandbox). Keep this window open (the display may
# turn off). The script refuses to start the pseudobulk with < 3 h to the 23:30 power cut, keeps the machine awake,
# and ends itself if its memory exceeds 4 GB. If it is interrupted, simply run it again: verified outputs are
# kept, anything incomplete is rebuilt.
param([string]$Sc = "H:\cs+bio_file\K562_gwps_raw_singlecell_01.h5ad")
$Root = Split-Path -Parent $PSScriptRoot
$EnvDir = Join-Path $env:USERPROFILE ".claude-science\conda\envs\pert-audit"
$env:PATH = "$EnvDir;$EnvDir\Library\bin;$EnvDir\Scripts;$env:PATH"
$env:PYTHONDONTWRITEBYTECODE = "1"
$Out = Join-Path $Root "results\conf_v0.7\gwps"
New-Item -ItemType Directory -Force $Out | Out-Null
$Log = Join-Path $Out ("stdout_c0_" + (Get-Date -Format "yyyyMMdd_HHmm") + ".txt")
Write-Host "GWPS C0 started $(Get-Date). Output: $Log"
& "$EnvDir\python.exe" (Join-Path $PSScriptRoot "confirm_run_v07.py") c0 --dataset gwps --sc $Sc *>&1 | Tee-Object -FilePath $Log
$rc = $LASTEXITCODE
Set-Content -Path (Join-Path $Out "c0_launcher_exit.json") -Value ("{""exit_code"": $rc, ""local_time"": ""$(Get-Date -Format s)""}")
Write-Host "confirm_run_v07.py c0 exited with code $rc (0 = done, 4 = QC failed, 10 = paused before the power cut, 7 = memory cap)"
