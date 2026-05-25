#!/usr/bin/env bash
# ===========================================================================
# EGAP_setup.sh
#
# Installation script for the Entheome Genome Assembly Pipeline (EGAP) v3.4.1.
# Installs Miniforge3 (if absent), creates the EGAP_env conda environment,
# installs auxiliary tools (runner, quast resources), and optionally provisions
# the Kraken2 16 GB Standard database used for read decontamination.
#
# Usage:
#   bash EGAP_setup.sh [--skip-kraken] [--kraken-prebuilt] [--threads N]
#                     [--kraken-db PATH]
#
# Options:
#   --skip-kraken      Skip the Kraken2 database step entirely.
#   --kraken-prebuilt  Download the pre-built 16 GB Standard database instead
#                      of building it from source. Faster but the filename may
#                      change over time — verify the latest at
#                      https://benlangmead.github.io/aws-indexes/k2
#   --threads N        Number of threads for Kraken2 build (default: 8).
#   --kraken-db PATH   Destination for the Kraken2 database
#                      (default: ${HOME}/kraken2_db).
#   -h, --help         Show this help text and exit.
#
# Environment variables set on success:
#   KRAKEN2_DB         Appended to ~/.bashrc and ~/.zshrc (if present)
#                      pointing at the chosen Kraken2 database directory.
#
# Author:  Ian Bollinger (ian.bollinger@entheome.org)
# Version: 3.4.1
# Updated: 2026-04-16
# ===========================================================================
set -euo pipefail


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
log_info()  { printf '\033[36m[INFO]\033[0m  %s\n'  "$*"; }
log_pass()  { printf '\033[32m[PASS]\033[0m  %s\n'  "$*"; }
log_warn()  { printf '\033[33m[WARN]\033[0m  %s\n'  "$*" >&2; }
log_error() { printf '\033[31m[ERROR]\033[0m %s\n'  "$*" >&2; }

die() {
    log_error "$*"
    exit 1
}


# ---------------------------------------------------------------------------
# Defaults and argument parsing
# ---------------------------------------------------------------------------
SKIP_KRAKEN=0
KRAKEN_PREBUILT=0
THREADS=8
KRAKEN2_DB_DIR="${HOME}/kraken2_db"
KRAKEN2_PREBUILT_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240904.tar.gz"

print_usage() {
    sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-kraken)      SKIP_KRAKEN=1; shift ;;
        --kraken-prebuilt)  KRAKEN_PREBUILT=1; shift ;;
        --threads)
            [[ $# -ge 2 ]] || die "--threads requires an integer argument"
            THREADS="$2"; shift 2 ;;
        --kraken-db)
            [[ $# -ge 2 ]] || die "--kraken-db requires a path argument"
            KRAKEN2_DB_DIR="$2"; shift 2 ;;
        -h|--help)          print_usage; exit 0 ;;
        *)                  die "Unknown option: $1 (use --help for usage)" ;;
    esac
done

[[ "${THREADS}" =~ ^[0-9]+$ ]] || die "--threads must be a positive integer (got '${THREADS}')"


# ---------------------------------------------------------------------------
# OS / architecture detection
# ---------------------------------------------------------------------------
detect_platform() {
    local uname_s uname_m
    uname_s="$(uname -s)"
    uname_m="$(uname -m)"

    case "${uname_s}" in
        Linux)  OS_NAME="Linux"  ;;
        Darwin) OS_NAME="MacOSX" ;;
        *)      die "Unsupported OS: ${uname_s}. EGAP supports Linux and macOS." ;;
    esac

    case "${uname_m}" in
        x86_64|amd64)         ARCH="x86_64" ;;
        aarch64|arm64)        ARCH="aarch64" ;;
        *)                    die "Unsupported architecture: ${uname_m}" ;;
    esac

    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${OS_NAME}-${ARCH}.sh"
    log_info "Detected platform: ${OS_NAME} / ${ARCH}"
}


# ---------------------------------------------------------------------------
# Miniforge3 installation
# ---------------------------------------------------------------------------
install_miniforge() {
    if command -v conda >/dev/null 2>&1; then
        log_pass "conda already on PATH: $(command -v conda)"
        return
    fi

    if [[ -d "${HOME}/miniforge3" ]]; then
        log_warn "Miniforge3 directory exists at ${HOME}/miniforge3 but conda is not on PATH."
        log_info "Sourcing existing installation."
    else
        log_info "Downloading Miniforge3 installer..."
        local installer="${HOME}/Miniforge3-installer.sh"
        wget -q "${MINIFORGE_URL}" -O "${installer}" \
            || die "Failed to download Miniforge3 installer from ${MINIFORGE_URL}"

        log_info "Installing Miniforge3 to ${HOME}/miniforge3..."
        bash "${installer}" -b -p "${HOME}/miniforge3" \
            || die "Miniforge3 installer returned a non-zero exit code"
        rm -f "${installer}"
    fi

    # shellcheck disable=SC1091
    source "${HOME}/miniforge3/etc/profile.d/conda.sh"
    "${HOME}/miniforge3/bin/conda" init bash >/dev/null 2>&1 || true
    log_pass "Miniforge3 ready."
}


# ---------------------------------------------------------------------------
# Create or update the EGAP_env conda environment
# ---------------------------------------------------------------------------
create_egap_env() {
    log_info "Creating/updating conda environment 'EGAP_env' (may take several minutes)..."

    if conda env list | awk '{print $1}' | grep -qx "EGAP_env"; then
        log_warn "EGAP_env already exists — updating the 'egap' package in place."
        conda install -y -n EGAP_env -c bioconda -c conda-forge "egap" \
            || die "Failed to update egap in existing EGAP_env"
    else
        conda create -y -n EGAP_env \
            "python>=3.8,<3.9" \
            egap \
            -c bioconda -c conda-forge \
            || die "Failed to create EGAP_env"
    fi

    log_pass "EGAP_env is ready."
}


# ---------------------------------------------------------------------------
# Install auxiliary tools inside EGAP_env
# ---------------------------------------------------------------------------
install_aux_tools() {
    log_info "Activating EGAP_env for auxiliary tool installation..."
    # shellcheck disable=SC1091
    source "${HOME}/miniforge3/etc/profile.d/conda.sh"
    conda activate EGAP_env || die "Could not activate EGAP_env"

    log_info "Cloning and installing 'runner' from github.com/dfguan/runner..."
    local tmpdir
    tmpdir="$(mktemp -d)"
    (
        cd "${tmpdir}"
        git clone --depth 1 https://github.com/dfguan/runner.git
        cd runner
        pip install --no-cache-dir .
    ) || die "runner installation failed"
    rm -rf "${tmpdir}"

    if python -c "import runner" >/dev/null 2>&1; then
        log_pass "runner installed successfully."
    else
        log_warn "runner import check failed — the package may not be functional."
    fi

    log_info "Downloading QUAST auxiliary datasets (gridss, silva)..."
    quast-download-gridss || log_warn "quast-download-gridss failed (non-fatal)"
    quast-download-silva  || log_warn "quast-download-silva failed (non-fatal)"
    log_pass "Auxiliary tools installed."
}


# ---------------------------------------------------------------------------
# Kraken2 database provisioning
# ---------------------------------------------------------------------------
kraken2_build_from_source() {
    log_info "Building Kraken2 standard 16 GB database from source at ${KRAKEN2_DB_DIR}"
    log_warn "This can take 6–12 hours and requires ~100 GB of free disk space."
    mkdir -p "${KRAKEN2_DB_DIR}"
    kraken2-build --standard --db "${KRAKEN2_DB_DIR}" --threads "${THREADS}" --use-ftp \
        || die "kraken2-build failed"
    log_pass "Kraken2 database built at ${KRAKEN2_DB_DIR}"
}

kraken2_download_prebuilt() {
    log_info "Downloading pre-built Kraken2 16 GB database to ${KRAKEN2_DB_DIR}"
    log_warn "Verify the URL is current at https://benlangmead.github.io/aws-indexes/k2"
    mkdir -p "${KRAKEN2_DB_DIR}"
    local tarball="${KRAKEN2_DB_DIR}/k2_standard_16gb.tar.gz"
    wget -c "${KRAKEN2_PREBUILT_URL}" -O "${tarball}" \
        || die "Failed to download pre-built Kraken2 database"
    log_info "Extracting Kraken2 archive (this may take 10-30 minutes)..."
    tar -xzf "${tarball}" -C "${KRAKEN2_DB_DIR}/" \
        || die "Failed to extract Kraken2 archive"
    rm -f "${tarball}"
    log_pass "Kraken2 database extracted to ${KRAKEN2_DB_DIR}"
}

provision_kraken2() {
    if [[ "${SKIP_KRAKEN}" -eq 1 ]]; then
        log_warn "Skipping Kraken2 setup (--skip-kraken). Read decontamination will be disabled until KRAKEN2_DB is set."
        return
    fi

    if [[ "${KRAKEN_PREBUILT}" -eq 1 ]]; then
        kraken2_download_prebuilt
    else
        kraken2_build_from_source
    fi

    export_kraken2_var
}

export_kraken2_var() {
    local export_line="export KRAKEN2_DB=\"${KRAKEN2_DB_DIR}\""
    for rc in "${HOME}/.bashrc" "${HOME}/.zshrc"; do
        [[ -f "${rc}" ]] || continue
        if grep -qF "${export_line}" "${rc}"; then
            log_info "KRAKEN2_DB already exported in ${rc}"
        else
            {
                echo ""
                echo "# Added by EGAP_setup.sh"
                echo "${export_line}"
            } >> "${rc}"
            log_pass "Appended KRAKEN2_DB export to ${rc}"
        fi
    done
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main() {
    log_info "EGAP v3.4.1 setup starting."
    detect_platform
    install_miniforge

    # shellcheck disable=SC1091
    source "${HOME}/miniforge3/etc/profile.d/conda.sh"

    create_egap_env
    install_aux_tools
    provision_kraken2

    log_pass "EGAP v3.4.1 installation complete."
    echo
    echo "Next steps:"
    echo "  1. Open a new shell (or run: source ~/.bashrc)"
    echo "  2. Activate the environment: conda activate EGAP_env"
    echo "  3. Run EGAP: EGAP --help"
    if [[ "${SKIP_KRAKEN}" -eq 0 ]]; then
        echo "  4. Kraken2 DB is at: ${KRAKEN2_DB_DIR}"
    else
        echo "  4. Set KRAKEN2_DB to your Kraken2 database path before running EGAP"
        echo "     (or decontamination will be skipped)."
    fi
}

main "$@"
