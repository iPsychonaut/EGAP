#!/bin/bash

# Function to check if a command was successful
check_success() {
    if [ $? -ne 0 ]; then
        echo -e "\e[31mERROR: $1 failed.\e[0m"
        exit 1
    fi
}

# Updating Java Runtime Environment (JRE)
echo -e "\e[36mUpdating Java Runtime Environment (JRE)...\e[0m"
sudo apt-get update && sudo apt install openjdk-8-jre-headless -y
check_success "JRE update"

# Installing unzip
echo -e "\e[36mDownloading unzip...\e[0m"
sudo apt install unzip
check_success "Unzip installation"

# Checking and Updating Conda
echo -e "\e[36mChecking for Conda...\e[0m"
if ! command -v conda &> /dev/null; then
    echo -e "\e[31mConda is not installed. Please install Conda before running this script.\e[0m"
    exit 1
else
    echo -e "\e[36mUpdating Conda...\e[0m"
    conda update -n base -c defaults conda -y
    check_success "Conda update"
fi

# Check if gdown is installed, and install it if it's not
if ! command -v gdown &> /dev/null; then
    echo -e "\e[36mInstalling gdown...\e[0m"
    pip install gdown==4.7.1
fi

# Check if file "EGAP_installs.zip" is in the current folder, if not download it from the following github link:
if [ ! -f "EGAP_installs.zip" ]; then
    echo -e "\e[36mDownloading EGAP_installs.zip...\e[0m"
    wget https://raw.githubusercontent.com/iPsychonaut/EGAP/126930051737949b2565e1daa9d0d40883fb0797/EGAP_installs.zip
    check_success "Downloading EGAP_installs.zip"
fi

# Unzip the contents of "EGAP_installs.zip" directly into this folder
echo -e "\e[36mUnzipping EGAP_installs.zip...\e[0m"
unzip -o EGAP_installs.zip
check_success "Unzipping EGAP_installs.zip"

# Remove EGAP_installs.zip
rm EGAP_installs.zip
check_success "Removing EGAP_installs.zip"

# Creating a custom Conda environment
echo -e "\e[36mCreating custom Conda environment 'entheome_env' with Python 3.8.15...\e[0m"
conda create -n entheome_env python=3.8.15 -y
check_success "Conda environment creation"

# Restart terminal, acivate entheome_env, and run mamba install command
echo -e "\e[32mInitial Setup complete.\e[0m"

echo -e "\e[36mRestart the terminal and run 'conda activate entheome_env' to activate the environment and then run the following 'mamba install...' command.\e[0m"

echo -e "\e[36mmamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel pandas==2.0.3 busco==5.5.0 nanoq==0.10.0 biopython==1.81 tqdm==4.38.0 beautifulsoup4==4.12.2 quast==5.2.0 nanostat==1.6.0 flye==2.9.2 bbtools==37.62 metaeuk==6.a5d39d9 blast==2.14.1 bwa==0.7.17 minimap2==2.26 pysam==0.21.0 samtools==1.17 arcs==1.2.5 tigmint==1.2.10 abyss==2.3.7 racon==1.5.0 spades==3.15.3 gdown==4.7.1 psutil==5.9.5 abyss==2.3.7 requests==2.31.0 minimap2==2.26 spoa==4.1.3 racon==1.5.0 termcolor==2.3.0 fastqc==0.12.1 masurca==4.1.0 openjdk=8\e[0m"
