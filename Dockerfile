# start with miniconda3 as build environment
FROM condaforge/mambaforge AS build
LABEL maintainer="Ian M Bollinger <ian.bollinger@entheome.org>"

# docker build -t entheome_ecosystem .
# docker run -it -v /mnt/d:/mnt/d entheome_ecosystem bash

# Update, install mamba and conda-pack:
RUN mamba install -n base --yes conda-pack

###############################################################################
# Generate EGAP_env
###############################################################################

# Install EGAP deps from bioconda
# here specifying specific versions to be able to set ENV below
RUN conda create -c conda-forge -c bioconda -c defaults \
    -n EGAP_env --yes python==3.8

RUN conda run -n EGAP_env mamba install -y -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    masurca==4.1.2 quast==5.2.0 compleasm==0.2.6 biopython==1.81 RagTag==2.1.0 \
    NanoPlot==1.43.0 termcolor==2.3.0 minimap2==2.28 bwa==0.7.18 matplotlib==3.7.3 \
    bamtools==2.5.2 tgsgapcloser==1.2.1 psutil==6.0.0 pilon==1.22 trimmomatic==0.39 \
    beautifulsoup4==4.12.3 ncbi-datasets-cli==16.39.0 samtools==1.21 fastqc==0.12.1 \
    abyss==2.0.2 bbmap==39.15 ratatosk==0.9.0 racon==1.5.0 kmc==3.2.4 runner==1.3 \
    spades==4.0.0 flye && \
    conda clean -a -y

# Download required resources for quast
RUN conda run -n EGAP_env quast-download-gridss && \
    conda run -n EGAP_env quast-download-silva

# Ensure make is available to build remaining libraries
RUN apt-get update && apt-get install -y make build-essential zlib1g-dev

# Clone purge_dups repository if it doesn't exist
RUN conda run -n EGAP_env git clone https://github.com/dfguan/purge_dups.git && \
                          cd purge_dups/src && \
    conda run -n EGAP_env make && \
                          cd ~

# Clone Ratatosk repository if it doesn't exist
RUN conda run -n EGAP_env git clone https://github.com/DecodeGenetics/Ratatosk.git && \
                          cd Ratatosk && \
                          mkdir build && \
                          cd build && \
    conda run -n EGAP_env cmake .. && \
    conda run -n EGAP_env make && \
    conda run -n EGAP_env make install && \
                          cd ~

# package with conda-pack
RUN conda-pack --ignore-missing-files -n EGAP_env -o /tmp/EGAP_env.tar && \
    mkdir /EGAP_env && cd /EGAP_env && tar xf /tmp/EGAP_env.tar && \
    rm /tmp/EGAP_env.tar && \
    /EGAP_env/bin/conda-unpack

###############################################################################
# Generate EGEP_env
###############################################################################

# Install EGEP deps from bioconda
# here specifying specific versions to be able to set ENV below
RUN conda create -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    -n EGEP_env --yes "python>=3.6,<3.9"


RUN conda run -n EGEP_env mamba install -y -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    ete3 distro glimmerhmm bamtools ucsc-pslcdnafilter trimmomatic raxml iqtree trimal hisat2 tantan bedtools tbl2asn blat hmmer exonerate minimap2 stringtie goatools matplotlib-base natsort numpy pigz pandas psutil requests scipy seaborn \
    perl perl-yaml perl-file-which perl-local-lib perl-dbd-mysql perl-clone perl-hash-merge \
    perl-soap-lite perl-json perl-logger-simple perl-scalar-util-numeric perl-math-utils perl-mce \
    perl-text-soundex perl-parallel-forkmanager perl-db-file perl-perl4-corelibs \
    "biopython<1.80" xlrd==1.2.0 "trinity==2.8.5" "evidencemodeler==1.1.1" "pasa==2.4.1" "codingquarry==2.0" "scikit-learn<1.0.0" "blast=2.2.31" "diamond>=2.0.5" "trnascan-se>=2.0" "mafft>=7" "kallisto==0.46.1" "salmon>=0.9" "samtools>=1.9" \
    proteinortho \
    && conda clean -a -y

RUN conda run -n EGEP_env mamba install -y -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    gffutils==0.12 \
    disjoint-set==0.7.4 \
    cblaster==1.3.18 \
    defusedxml==0.7.1 \
    orthofinder==3.0.1b1 && \
    conda clean -a -y
    
# get necessary pip installs
SHELL ["conda", "run", "-n", "EGEP_env", "/bin/bash", "-c"]
RUN pip install clinker genomicsqlite

# Get the most recent version of Funannotate
SHELL ["conda", "run", "-n", "EGEP_env", "/bin/bash", "-c"]
RUN python -m pip install git+https://github.com/nextgenusfs/funannotate.git

# package with conda-pack
RUN conda-pack --ignore-missing-files -n EGEP_env -o /tmp/EGEP_env.tar && \
    mkdir /EGEP_env && cd /EGEP_env && tar xf /tmp/EGEP_env.tar && \
    rm /tmp/EGEP_env.tar && \
    /EGEP_env/bin/conda-unpack

###############################################################################
# Create Final Debian Build for Augustus
###############################################################################

# Now build environment
FROM debian:buster AS runtime

# Copy BOTH conda envs from the build stage
COPY --from=build /EGEP_env /EGEP_env
COPY --from=build /EGAP_env /EGAP_env

# Install Debian packages
RUN apt-get update && apt-get install -y \
    snap augustus augustus-data locales locales-all libgl1 procps \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/snap-hmm /usr/bin/snap \
    && rm -f "/EGEP_env/bin/fasta" \
    && ln -s "/EGEP_env/bin/fasta36" "/EGEP_env/bin/fasta"

# Add either or both to PATH
ENV PATH="/EGEP_env/bin:/EGAP_env/bin:$PATH"
ENV AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config" \
    EVM_HOME="/EGEP_env/opt/evidencemodeler-1.1.1" \
    PASAHOME="/EGEP_env/opt/pasa-2.4.1" \
    TRINITYHOME="/EGEP_env/opt/trinity-2.8.5" \
    QUARRY_PATH="/EGEP_env/opt/codingquarry-2.0/QuarryFiles" \
    ZOE="/usr/share/snap" \
    USER="me" \
    FUNANNOTATE_DB="/opt/databases" \
    GENEMARK_PATH="/mnt/d/EGEP" # REQUIRES USER TO DOWNLOAD AND UPDATE PATH

RUN conda run -n EGEP_env funannotate setup -d "/opt/databases"

# SHELL ["/bin/bash", "-c"]
# source /EGAP_env/bin/activate
# source /EGEP_env/bin/activate
