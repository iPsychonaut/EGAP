# Make the external E: drive (exFAT, holds E:\egap-ab) visible to the k3d cluster.
#
# Docker Desktop's WSL2 backend only auto-mounts fixed NTFS drives into its VM.
# Removable/exFAT drives get a placeholder /mnt/host/e with nothing mounted, so a
# bind of E:\egap-ab silently becomes an empty scratch directory and the k3d
# node's /data is empty. The VM-side mount does not survive a Docker Desktop
# restart, so run this after Docker starts (and after plugging the drive in),
# then restart the node so its /data bind re-resolves.
#
#   powershell -File tests\ab\k3d-mount-e.ps1
#
# ponytail: a manual step, not a Docker/WSL startup hook. Automate it if the
# drive stays attached permanently and the step keeps getting forgotten.

$ErrorActionPreference = 'Stop'

if (-not (Test-Path 'E:\egap-ab')) { throw 'E:\egap-ab not found. Is the external drive plugged in?' }

wsl -d docker-desktop -e sh -c 'mountpoint -q /mnt/host/e || mount -t drvfs E: /mnt/host/e'
if ($LASTEXITCODE -ne 0) { throw 'mounting E: inside the docker-desktop VM failed' }

docker restart k3d-egap-server-0 | Out-Null

# Prove it: the node must see the A/B inputs.
$seen = docker exec k3d-egap-server-0 sh -c 'ls /data/input 2>/dev/null | head -1'
if (-not $seen) { throw 'k3d node /data is still empty after the remount' }
Write-Host "OK: k3d node /data -> E:\egap-ab (sees $seen)"
