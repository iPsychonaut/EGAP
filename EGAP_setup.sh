#!/bin/bash

# Change to the home directory
cd ~ || { echo "Failed to change to home directory."; exit 1; }

# Function to check if a command was successful
check_success() {
                 if [ $? -ne 0 ]; then
                    echo -e "\e[31mERROR: $1 failed.\e[0m"
                    exit 1
                 fi
                }

# Download and install Miniforge3 if not already installed
if [ -d "$HOME/miniforge3" ]; then
    echo -e "\e[33mMiniforge3 is already installed at $HOME/miniforge3. Skipping installation.\e[0m"
else
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    echo -e "\e[36m\nDownloading Miniforge3 installer...\e[0m"
    wget "$MINIFORGE_URL" -O Miniforge3-Linux-x86_64.sh
    check_success "Downloading Miniforge3-Linux-x86_64.sh"

    # Run install bash script
    echo -e "\e[36m\nInstalling Miniforge3...\e[0m"
    bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
    check_success "Installing Miniforge3"

    # Cleanup Miniforge3 installer
    rm Miniforge3-Linux-x86_64.sh
    check_success "Removing Miniforge3 installer"

    # Initialize conda
    "$HOME/miniforge3/bin/conda" init
    check_success "Initializing conda"
fi

# Initialize mamba for the current shell
echo -e "\e[36m\nInitializing mamba for the current shell...\e[0m"
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
check_success "Initializing mamba"

# Source .bashrc to update shell environment
source ~/.bashrc
check_success "\nSourcing ~/.bashrc"

# Create EGAP_env with Python 3.8 using mamba
echo -e "\e[36m\nCreating conda environment 'EGAP_env' with Python 3.8...\e[0m"
conda create -y -n EGAP_env python=3.8
check_success "Creating EGAP_env"

# Activate the environment
echo -e "\e[36m\nActivating 'EGAP_env' environment...\e[0m"
conda activate EGAP_env
check_success "Activating 'EGAP_env' environment"

# EGAP Installs
echo -e "\e[36m\nInstalling system packages via apt-get...\e[0m"
sudo apt-get update && sudo apt-get install -y openjdk-8-jre-headless
check_success "Installing system packages"

# Install required conda packages
echo -e "\e[36m\nInstalling required conda packages via mamba...\e[0m"
mamba install -y -c bioconda -c conda-forge \
                                masurca==4.1.2 \
                                quast==5.2.0 \
                                compleasm==0.2.6 \
                                biopython==1.81 \
                                RagTag==2.1.0 \
                                NanoPlot==1.43.0 \
                                termcolor==2.3.0 \
                                minimap2==2.28 \
                                bwa==0.7.18 \
                                samtools==1.21 \
                                bamtools==2.5.2 \
                                tgsgapcloser==1.2.1 \
                                abyss==2.0.2 \
                                sepp==4.5.1 \
                                psutil==6.0.0 \
                                beautifulsoup4==4.12.3 \
                                ncbi-datasets-cli==16.39.0 \
                                matplotlib==3.7.3 \
                                trimmomatic==0.39 \
                                pilon==1.22 \
                                fastqc==0.12.1 \
                                bbmap==39.15 \
                                racon==1.5.0 \
                                kmc==3.2.4 \
                                runner==1.3 \
                                spades==4.0.0 \
                                flye

check_success "Installing conda packages"

# Clone purge_dups repository if it doesn't exist
if [ ! -d "purge_dups" ]; then
    echo -e "\e[36m\nCloning purge_dups repository...\e[0m"
    git clone https://github.com/dfguan/purge_dups.git
    check_success "Cloning purge_dups repository"
else
    echo -e "\e[33mDirectory 'purge_dups' already exists. Skipping clone.\e[0m"
fi

# Build purge_dups
cd purge_dups/src || { echo "Failed to change to purge_dups/src directory."; exit 1; }
make
check_success "Building purge_dups"
cd ~ # Return to home directory

# Clone Ratatosk repository if it doesn't exist
if [ ! -d "Ratatosk" ]; then
    echo -e "\e[36m\nCloning Ratatosk repository...\e[0m"
    git clone --recursive https://github.com/DecodeGenetics/Ratatosk.git
    check_success "Cloning Ratatosk repository"
else
    echo -e "\e[33mDirectory 'Ratatosk' already exists. Skipping clone.\e[0m"
fi

# Build Ratatosk
cd Ratatosk || { echo "Failed to change to Ratatosk directory."; exit 1; }
mkdir -p build && cd build
cmake .. || { echo "CMake configuration failed."; exit 1; }
make || { echo "Make failed."; exit 1; }
sudo make install || { echo "Make install failed."; exit 1; }
check_success "Building Ratatosk"
cd ~ # Return to home directory

echo -e "\n\e[32mEGAP Pipeline Pre-requisites have been successfully installed!\e[0m\n"
echo -e '\n\e[32mStart by activating the environment "conda activate EGAP_env"\e[0m\n'
echo -e '\n\e[36mOptionally run "quast-download-gridss"\e[0m\n'
echo -e '\n\e[36mOptionally run "quast-download-silva"\e[0m\n'
