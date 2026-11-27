#!/usr/bin/env bash
set -eux

# Install the same tree the git checkout and the Docker image use:
#   $PREFIX/share/egap/EGAP.py
#   $PREFIX/share/egap/bin/*.py
#   $PREFIX/share/egap/resources/...
# EGAP.py finds bin/ beside itself, and bin/ finds ../resources, so no layout
# needs its own lookup and nothing is downloaded at run time.
SHARE="${PREFIX}/share/egap"
mkdir -p "${SHARE}/bin" "${SHARE}/resources" "${PREFIX}/bin"
cp EGAP.py "${SHARE}/EGAP.py"
cp bin/*.py "${SHARE}/bin/"
cp -r resources/templates "${SHARE}/resources/templates"
cp resources/TruSeq3-PE.fa "${SHARE}/resources/TruSeq3-PE.fa"

# Launchers. Python rather than shell so the noarch package runs anywhere, and
# located relative to the launcher itself so they work from any prefix.
for entry in "EGAP:EGAP.py" "EGAP_TUI:bin/EGAP_TUI.py"; do
    command_name="${entry%%:*}"
    script_path="${entry##*:}"
    cat > "${PREFIX}/bin/${command_name}" <<EOF
#!/usr/bin/env python3
import os
import runpy

prefix = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
runpy.run_path(os.path.join(prefix, "share", "egap", "${script_path}"), run_name="__main__")
EOF
    chmod +x "${PREFIX}/bin/${command_name}"
done
