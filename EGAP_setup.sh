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

# Check if file "EGAP_installs.zip" is in the current folder, if not download it from the following GitHub link:
if [ ! -f "EGAP_installs.zip" ]; then
    echo -e "\e[36mDownloading EGAP_installs.zip...\e[0m"
    wget https://raw.githubusercontent.com/iPsychonaut/EGAP/126930051737949b2565e1daa9d0d40883fb0797/EGAP_installs.zip
    check_success "Downloading EGAP_installs.zip"
fi

# Add downloaded tools to PATH if not already present
if ! grep -q 'export PATH="$HOME/Pilon-1.24:$PATH"' ~/.bashrc; then
    echo 'export PATH="$HOME/Pilon-1.24:$PATH"' >> ~/.bashrc
fi

if ! grep -q 'export PATH="$HOME/Trimmomatic-0.39:$PATH"' ~/.bashrc; then
    echo 'export PATH="$HOME/Trimmomatic-0.39:$PATH"' >> ~/.bashrc
fi

# Unzip the contents of "EGAP_installs.zip" directly into this folder
echo -e "\e[36mUnzipping EGAP_installs.zip...\e[0m"
unzip -o EGAP_installs.zip
check_success "Unzipping EGAP_installs.zip"

# Remove EGAP_installs.zip
rm EGAP_installs.zip
check_success "Removing EGAP_installs.zip"

# Clone purge_dups repository if it doesn't exist
if [ ! -d "purge_dups" ]; then
    echo -e "\e[36mCloning purge_dups repository...\e[0m"
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

# Download and install Miniforge3 if not already installed
if [ -d "$HOME/miniforge3" ]; then
    echo -e "\e[33mMiniforge3 is already installed at $HOME/miniforge3. Skipping installation.\e[0m"
    # Initialize Conda/Mamba
    eval "$(conda shell.bash hook)"
else
    MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    echo -e "\e[36mDownloading Miniforge3 installer...\e[0m"
    wget "$MINIFORGE_URL" -O Miniforge3-Linux-x86_64.sh
    check_success "Downloading Miniforge3-Linux-x86_64.sh"

    # Run install bash script
    echo -e "\e[36mInstalling Miniforge3...\e[0m"
    bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
    check_success "Installing Miniforge3"

    # Cleanup Miniforge3 installer
    rm Miniforge3-Linux-x86_64.sh
    check_success "Removing Miniforge3 installer"
    
    # Initialize conda
    "$HOME/miniforge3/bin/conda" init
    check_success "Initializing conda"
    
    # Initialize Conda/Mamba
    eval "$(conda shell.bash hook)"
    
    # Source .bashrc to update PATH
    source ~/.bashrc
    check_success "Sourcing .bashrc"
fi

# Install mamba into the base environment
echo -e "\e[36mInstalling mamba into the base Conda environment...\e[0m"
conda install -y -c conda-forge mamba
check_success "Installing mamba"

# Create entheome_env with Python 3.8 using mamba
echo -e "\e[36mCreating conda environment 'entheome_env' with Python 3.8...\e[0m"
mamba create -y -n entheome_env python=3.8
check_success "Creating entheome_env"

# EGAP Installs
echo -e "\e[36mInstalling system packages via apt-get...\e[0m"
sudo apt-get update && sudo apt-get install -y openjdk-8-jre-headless racon fastqc
check_success "Installing system packages"

# Install required conda packages
echo -e "\e[36mInstalling required conda packages via mamba...\e[0m"
conda activate entheome_env && mamba install -y -c bioconda -c conda-forge \
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
                                                               abyss==2.3.10 \
                                                               sepp==4.5.1 \
                                                               psutil==6.0.0 \
                                                               merqury==1.3 \
                                                               meryl==1.3 \
                                                               beautifulsoup4

check_success "Installing conda packages"

# # TODO: Get QUAST follow-up commands to function appropriately
# # QUAST Support Downloads
# echo -e "\e[36mDownloading QUAST support data...\e[0m"
# # Ensure QUAST's bin directory is in PATH
# export PATH="$HOME/miniforge3/envs/entheome_env/bin:$PATH"
# quast-download-grids
# check_success "Downloading QUAST grids"
# quast-download-silva
# check_success "Downloading QUAST silva"

echo -e "\n\e[32mEGAP Pipeline Pre-requisites have been successfully installed!\e[0m\n"
echo -e '\n\e[32mstart by activating the environment "mamba activate entheome_env"\e[0m\n'
echo -e '\n\e[36Optionally run "quast-download-grids"\e[0m\n'
echo -e '\n\e[36Optionally run "quast-download-silva"\e[0m\n'
