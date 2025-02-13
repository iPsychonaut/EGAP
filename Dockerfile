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

# Create EGAP_env with Python 3.8
# RUN conda create -c conda-forge -c bioconda -c defaults -n EGAP_env --yes python==3.8 egap
RUN conda create -c conda-forge -c bioconda -c defaults \
    -n EGAP_env --yes python=3.8 egap
    
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

# Download the latest EGAP.py from GitHub into the EGAP_env.
# Install wget (if not already available) to retrieve the file.
RUN apt-get update && apt-get install -y wget && \
    wget -O /EGAP_env/EGAP.py https://raw.githubusercontent.com/iPsychonaut/EGAP/master/EGAP.py && \
    chmod +x /EGAP_env/EGAP.py && \
    rm -rf /var/lib/apt/lists/*

# Create a wrapper script in EGAP_env/bin called "EGAP"
RUN printf '#!/bin/bash\npython /EGAP_env/EGAP.py "$@"\n' > /EGAP_env/bin/EGAP && \
    chmod +x /EGAP_env/bin/EGAP

###############################################################################
# Generate EGEP_env
###############################################################################

# Create EGEP_env with Python (version between 3.6 and 3.9)
RUN conda create -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    -n EGEP_env --yes "python>=3.6,<3.9"

# Install EGEP dependencies from bioconda
RUN conda run -n EGEP_env mamba install -y -c conda-forge -c bioconda -c prkrekel -c agbiome -c defaults \
    perl perl-yaml perl-file-which perl-local-lib perl-dbd-mysql perl-clone perl-hash-merge \
    perl-soap-lite perl-json perl-logger-simple perl-scalar-util-numeric perl-math-utils perl-mce \
    perl-text-soundex perl-parallel-forkmanager perl-db-file perl-perl4-corelibs \
    ete3 distro glimmerhmm bamtools ucsc-pslcdnafilter trimmomatic raxml iqtree trimal hisat2 tantan bedtools tbl2asn blat hmmer exonerate minimap2 stringtie goatools matplotlib-base natsort numpy pigz pandas psutil requests scipy seaborn proteinortho

RUN conda run -n EGEP_env mamba install -y -c conda-forge -c bioconda -c prkrekel -c defaults \
    gffutils==0.12 disjoint-set==0.7.4 cblaster==1.3.18  defusedxml==0.7.1 orthofinder==3.0.1b1 \
    "biopython<1.80" xlrd==1.2.0 "trinity==2.8.5" "evidencemodeler==1.1.1" "pasa==2.4.1" "codingquarry==2.0" "scikit-learn<1.0.0" "blast=2.2.31" "diamond>=2.0.5" "trnascan-se>=2.0" "mafft>=7" "kallisto==0.46.1" "salmon>=0.9" "samtools>=1.9" \
    && conda clean -a -y
    
# Get necessary pip installs for EGEP_env
SHELL ["conda", "run", "-n", "EGEP_env", "/bin/bash", "-c"]
RUN pip install clinker genomicsqlite

# Get the most recent version of Funannotate
SHELL ["conda", "run", "-n", "EGEP_env", "/bin/bash", "-c"]
RUN python -m pip install git+https://github.com/nextgenusfs/funannotate.git

# Package EGEP_env with conda-pack
RUN conda-pack --ignore-missing-files -n EGEP_env -o /tmp/EGEP_env.tar && \
    mkdir /EGEP_env && cd /EGEP_env && tar xf /tmp/EGEP_env.tar && \
    rm /tmp/EGEP_env.tar && \
    /EGEP_env/bin/conda-unpack

###############################################################################
# Create Final Debian Build for Augustus
###############################################################################

# Build runtime image
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

# Set PATH so both environments' bin directories are available
ENV PATH="/EGEP_env/bin:/EGAP_env/bin:$PATH"
ENV AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config" \
    EVM_HOME="/EGEP_env/opt/evidencemodeler-1.1.1" \
    PASAHOME="/EGEP_env/opt/pasa-2.4.1" \
    TRINITYHOME="/EGEP_env/opt/trinity-2.8.5" \
    QUARRY_PATH="/EGEP_env/opt/codingquarry-2.0/QuarryFiles" \
    ZOE="/usr/share/snap" \
    USER="me" \
    FUNANNOTATE_DB="/opt/databases" \
    GENEMARK_PATH="/mnt/d/EGEP"

# Run funannotate setup (within EGEP_env)
RUN /EGEP_env/bin/funannotate setup -d "/opt/databases"
