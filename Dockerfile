# =============================================================================
# Entheome Ecosystem — multi-stage Docker image
# -----------------------------------------------------------------------------
# Builds three conda environments used across the Entheome toolchain:
#   * EGAP_env         — Entheome Genome Assembly Pipeline (EGAP) v3.4.0
#   * EGEP_env         — Annotation helper tools
#   * funannotate_env  — Funannotate for eukaryotic genome annotation
# The final runtime image is a slim Debian layer that carries all three envs
# plus Augustus and related runtime dependencies.
# =============================================================================

FROM condaforge/mambaforge AS build

LABEL maintainer="Ian Bollinger <ian.bollinger@entheome.org>" \
      version="3.4.0" \
      description="Entheome Genome Assembly Pipeline (EGAP) v3.4.0 — multi-env Entheome ecosystem image" \
      org.opencontainers.image.source="https://github.com/iPsychonaut/EGAP" \
      org.opencontainers.image.version="3.4.0" \
      org.opencontainers.image.authors="Ian Bollinger <ian.bollinger@entheome.org>"

# Install mamba and conda-pack in the base env — used to build and then
# pack per-env tarballs that are copied into the slim runtime image below.
RUN mamba install -n base --yes conda-pack

###############################################################################
# Generate EGAP_env
###############################################################################

# Create EGAP_env with Python 3.8 and all EGAP v3.4.0 dependencies.
# Pinning rationale (verified against conda list on 2026-04-06):
#   numpy=1.19.5     — tiara=1.0.3 requires numpy<1.20; do not loosen.
#   tiara=1.0.3      — exact pin; newer solves break the numpy constraint.
#   kraken2=2.1.6    — exact pin; tested working version.
#   sra-tools=3.2.0  — exact pin; API changes between minor versions.
#   trimmomatic=0.40 — exact pin; share-dir path used in adapter symlinks.
#   flye=2.9.5       — exact pin; >=3.0 changes assembly graph format.
RUN conda create -n EGAP_env -y -c bioconda -c conda-forge \
    'python>=3.8,<3.9' \
    pandas \
    'numpy=1.19.5' \
    'masurca=4.1.4' \
    'quast=5.3.0' \
    compleasm \
    busco \
    biopython \
    ragtag \
    'nanoplot=1.46.2' \
    termcolor \
    minimap2 \
    bwa-mem2 \
    samtools \
    bamtools \
    tgsgapcloser \
    abyss \
    sepp \
    psutil \
    beautifulsoup4 \
    ncbi-datasets-cli \
    matplotlib-base \
    'trimmomatic=0.40' \
    pilon \
    fastqc \
    bbmap \
    racon \
    kmc \
    'spades=4.2.0' \
    purge_dups \
    'flye=2.9.5' \
    pbccs \
    'hifiasm=0.25.0' \
    gfatools \
    bifrost \
    ratatosk \
    'sra-tools=3.2.0' \
    filtlong \
    pyinaturalist \
    jinja2 \
    geopy \
    tabulate \
    openpyxl \
    requests \
    rich \
    textual \
    'tiara=1.0.3' \
    'kraken2=2.1.6'
    
# Download required resources for quast
RUN conda run -n EGAP_env quast-download-gridss && \
    conda run -n EGAP_env quast-download-silva

# Clone and install runner in EGAP_env
RUN git clone https://github.com/dfguan/runner.git && \
    cd runner && conda run -n EGAP_env python3 setup.py install --user && \
    cd .. && rm -rf runner

# Package EGAP_env with conda-pack
RUN conda-pack --ignore-missing-files -n EGAP_env -o /tmp/EGAP_env.tar && \
    mkdir /EGAP_env && cd /EGAP_env && tar xf /tmp/EGAP_env.tar && \
    rm /tmp/EGAP_env.tar && \
    /EGAP_env/bin/conda-unpack

# Download EGAP v3.4.0 scripts from GitHub into the EGAP_env.
# Install wget (if not already available) to retrieve the files.
RUN apt-get update && apt-get install -y wget && \
    EGAP_BRANCH="v3.4.0" && \
    EGAP_RAW="https://raw.githubusercontent.com/iPsychonaut/EGAP/${EGAP_BRANCH}" && \
    wget -O /EGAP_env/EGAP.py "${EGAP_RAW}/EGAP.py" && \
    chmod +x /EGAP_env/EGAP.py && \
    mkdir -p /EGAP_env/bin && \
    for SCRIPT in \
        EGAP_TUI.py \
        file_manager.py \
        utilities.py \
        preprocess_illumina.py \
        preprocess_ont.py \
        preprocess_pacbio.py \
        preprocess_refseq.py \
        decontaminate_reads.py \
        decontaminate_assembly.py \
        assemble_masurca.py \
        assemble_spades.py \
        assemble_flye.py \
        assemble_hifiasm.py \
        compare_assemblies.py \
        polish_assembly.py \
        curate_assembly.py \
        qc_assessment.py \
        html_reporter.py \
        process_metadata.py \
        final_compress.py; \
    do \
        wget -O "/EGAP_env/bin/${SCRIPT}" "${EGAP_RAW}/bin/${SCRIPT}"; \
    done && \
    chmod +x /EGAP_env/bin/*.py && \
    rm -rf /var/lib/apt/lists/*

# Create a wrapper script in EGAP_env/bin called "EGAP"
RUN printf '#!/bin/bash\npython /EGAP_env/EGAP.py "$@"\n' > /EGAP_env/bin/EGAP && \
    chmod +x /EGAP_env/bin/EGAP

###############################################################################
# Generate EGEP_env
###############################################################################

RUN apt-get update && apt-get install -y \
    wget curl bzip2 ncbi-blast+ git build-essential autoconf automake libtool \
    pkg-config libglib2.0-dev zlib1g-dev snap augustus augustus-data locales \
    locales-all libgl1 procps && \
    rm -rf /var/lib/apt/lists/*

# Create EGEP_env with Python (version greater than 3.10)
RUN conda create -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    -n EGEP_env --yes python>=3.10

# Install EGEP dependencies from bioconda
RUN conda run -n EGEP_env mamba install -y -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    perl perl-yaml perl-file-which perl-local-lib perl-dbd-mysql perl-clone perl-hash-merge \
    perl-soap-lite perl-json perl-logger-simple perl-scalar-util-numeric perl-math-utils perl-mce \
    perl-text-soundex perl-parallel-forkmanager perl-db-file perl-perl4-corelibs \
    ete3 distro glimmerhmm bamtools ucsc-pslcdnafilter trimmomatic raxml iqtree trimal \
    hisat2 tantan bedtools tbl2asn blat hmmer exonerate minimap2 stringtie goatools \
    matplotlib-base natsort numpy pigz pandas psutil requests scipy seaborn proteinortho \
    clinker-py gffutils disjoint-set cblaster defusedxml orthofinder \
    biopython xlrd trinity evidencemodeler pasa codingquarry scikit-learn blast \
    diamond trnascan-se mafft kallisto salmon samtools \
    && conda clean -a -y
    
# Get necessary pip installs for EGEP_env
SHELL ["conda", "run", "-n", "EGEP_env", "/bin/bash", "-c"]
RUN pip install genomicsqlite

# Package EGEP_env with conda-pack
RUN conda-pack --ignore-missing-files -n EGEP_env -o /tmp/EGEP_env.tar && \
    mkdir /EGEP_env && cd /EGEP_env && tar xf /tmp/EGEP_env.tar && \
    rm /tmp/EGEP_env.tar && \
    /EGEP_env/bin/conda-unpack

###############################################################################
# Generate funannotate_env
###############################################################################

# Create funannotate_env with Python (version between 3.6 and 3.9)
RUN conda create -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    -n funannotate_env --yes "python>=3.6,<3.9"

RUN conda run -n funannotate_env mamba install -y -c conda-forge -c bioconda -c prkrekel -c defaults \
    gffutils==0.12 disjoint-set==0.7.4 cblaster==1.3.18  defusedxml==0.7.1 orthofinder==3.0.1b1 \
    "biopython<1.80" xlrd==1.2.0 "trinity==2.8.5" "evidencemodeler==1.1.1" "pasa==2.4.1" \
    "codingquarry==2.0" "scikit-learn<1.0.0" "blast=2.2.31" "diamond>=2.0.5" "trnascan-se>=2.0" \
    "mafft>=7" "kallisto==0.46.1" "salmon>=0.9" "samtools>=1.9" \
    && conda clean -a -y

# Get the most recent version of Funannotate
SHELL ["conda", "run", "-n", "funannotate_env", "/bin/bash", "-c"]
RUN python -m pip install git+https://github.com/nextgenusfs/funannotate.git

# Package funannotate_env with conda-pack
RUN conda-pack --ignore-missing-files -n EGEP_env -o /tmp/funannotate_env.tar && \
    mkdir /funannotate_env && cd /funannotate_env && tar xf /tmp/funannotate_env.tar && \
    rm /tmp/funannotate_env.tar && \
    /funannotate_env/bin/conda-unpack

###############################################################################
# Create Final Debian Build for Augustus
###############################################################################

# Build runtime image
FROM debian:buster AS runtime

# Copy BOTH conda envs from the build stage
COPY --from=build /EGAP_env /EGAP_env
COPY --from=build /EGEP_env /EGEP_env
COPY --from=build /funannotate_env /funannotate_env

# Install Debian packages
RUN apt-get update && apt-get install -y \
    snap augustus augustus-data locales locales-all libgl1 procps \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/snap-hmm /usr/bin/snap \
    && rm -f "/funannotate_env/bin/fasta" \
    && ln -s "/funannotate_env/bin/fasta36" "/funannotate_env/bin/fasta"

# Set PATH so both environments' bin directories are available
ENV PATH="/EGEP_env/bin:/EGAP_env/bin:/funannotate_env/bin:$PATH"
ENV AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config" \
    EVM_HOME="/funannotate_env/opt/evidencemodeler-1.1.1" \
    PASAHOME="/funannotate_env/opt/pasa-2.4.1" \
    TRINITYHOME="/funannotate_env/opt/trinity-2.8.5" \
    QUARRY_PATH="/funannotate_env/opt/codingquarry-2.0/QuarryFiles" \
    ZOE="/usr/share/snap" \
    USER="me" \
    FUNANNOTATE_DB="/opt/databases" \
    GENEMARK_PATH="/mnt/d/EGEP"

# -----------------------------------------------------------------------------
# EGAP v3.4.0 runtime defaults
# -----------------------------------------------------------------------------
# Kraken2 DB is NOT baked into the image (the standard 16 GB archive would
# roughly double the image size). Bind-mount the database directory at runtime
# and point KRAKEN2_DB at the mount point, e.g.:
#
#   docker run --rm \
#       -e KRAKEN2_DB=/kraken2_db \
#       -v /host/path/to/kraken2_db:/kraken2_db:ro \
#       -v /host/data:/data \
#       entheome_ecosystem:3.4.0 --input /data/samples.csv --output /data/out
#
# To provision a Kraken2 database on the host before running EGAP, either
# build from source (authoritative, ~6-12 hrs):
#   kraken2-build --standard --db /host/kraken2_db --threads 16 --use-ftp
# or download the pre-built standard 16 GB index (faster, verify filename at
# https://benlangmead.github.io/aws-indexes/k2):
#   wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240904.tar.gz \
#        -O /host/kraken2_db/k2_standard_16gb.tar.gz
#   tar -xzf /host/kraken2_db/k2_standard_16gb.tar.gz -C /host/kraken2_db/
ENV KRAKEN2_DB="" \
    CONDA_DEFAULT_ENV=EGAP_env \
    PYTHONUNBUFFERED=1 \
    EGAP_DRY_RUN=0

# One-time Funannotate database setup. Skipped silently if already initialised.
RUN /funannotate_env/bin/funannotate setup -d "/opt/databases"

# -----------------------------------------------------------------------------
# Default entrypoint — runs EGAP directly.
# Override to enter an interactive shell:
#   docker run --rm -it --entrypoint bash entheome_ecosystem:3.4.0
# -----------------------------------------------------------------------------
ENTRYPOINT ["/EGAP_env/bin/EGAP"]
CMD ["--help"]

# =============================================================================
# Usage examples
# -----------------------------------------------------------------------------
# Build:
#   docker build -t entheome_ecosystem:3.4.0 .
#
# Show EGAP help:
#   docker run --rm entheome_ecosystem:3.4.0
#
# Run EGAP with a host-mounted Kraken2 DB and data directory:
#   docker run --rm \
#       -e KRAKEN2_DB=/kraken2_db \
#       -v /host/kraken2_db:/kraken2_db:ro \
#       -v /host/data:/data \
#       entheome_ecosystem:3.4.0 \
#       --input-csv /data/samples.csv \
#       --output /data/output \
#       --threads 16 --ram 64
#
# Interactive shell (all three conda envs available on PATH):
#   docker run --rm -it --entrypoint bash entheome_ecosystem:3.4.0
# =============================================================================
