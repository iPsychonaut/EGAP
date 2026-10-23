# Make the external E: drive (exFAT, holds E:\egap-ab) visible to the k3d cluster.
#
# Docker Desktop's WSL2 VM mounts E: at /mnt/host/e only if the drive is attached
# when Docker starts (verified 2026-10-23: a Docker restart with the drive in
# brings it back by itself). Plug the drive in after Docker is up and you get an
# empty placeholder instead: a bind of E:\egap-ab silently becomes an empty
# scratch directory and the k3d node's /data is empty. Run this in that case
# (or whenever a pod reports an empty /data); it mounts the drive in the VM and
# restarts the node so its /data bind re-resolves.
#
#   powershell -File tests\ab\k3d-mount-e.ps1
#
# ponytail: a manual step, not a Docker/WSL startup hook. Automate it if the
# drive is hot-plugged often enough that the step keeps getting forgotten.

$ErrorActionPreference = 'Stop'

if (-not (Test-Path 'E:\egap-ab')) { throw 'E:\egap-ab not found. Is the external drive plugged in?' }

wsl -d docker-desktop -e sh -c 'mountpoint -q /mnt/host/e || mount -t drvfs E: /mnt/host/e'
if ($LASTEXITCODE -ne 0) { throw 'mounting E: inside the docker-desktop VM failed' }

docker restart k3d-egap-server-0 | Out-Null

# Prove it: the node must see the A/B inputs.
$seen = docker exec k3d-egap-server-0 sh -c 'ls /data/input 2>/dev/null | head -1'
if (-not $seen) { throw 'k3d node /data is still empty after the remount' }
Write-Host "OK: k3d node /data -> E:\egap-ab (sees $seen)"
