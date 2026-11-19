# run_gwps_download.ps1 -- user-run launcher for the K562 GWPS raw single-cell download.
# Runs outside the sandbox (H: is exFAT, the sandbox cannot write there). Keep this window open;
# the display may turn off. To pause cleanly: create the file <Dest>.stop (or press Ctrl+C; the
# download resumes from the state file next time).
param(
    [string]$Dest = "H:\cs+bio_file\K562_gwps_raw_singlecell_01.h5ad",
    [int]$Threads = 16,
    [double]$ChunkMiB = 32
)
$Root = Split-Path -Parent $PSScriptRoot                    # project root (parent of code\)
$EnvDir = Join-Path $env:USERPROFILE ".claude-science\conda\envs\pert-audit"
$env:PATH = "$EnvDir;$EnvDir\Library\bin;$EnvDir\Scripts;$env:PATH"
$env:PYTHONDONTWRITEBYTECODE = "1"
$Status = Join-Path $Root "logs\gwps_download_status.json"
& "$EnvDir\python.exe" (Join-Path $PSScriptRoot "download_resumable.py") `
    --file-id 35775507 --dest $Dest --threads $Threads --chunk-mib $ChunkMiB --status-file $Status
Write-Host "download_resumable.py exited with code $LASTEXITCODE"
