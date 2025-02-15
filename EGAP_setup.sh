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

# Create EGAP_env with Python using mamba
echo -e "\e[36m\nCreating conda environment 'EGAP_env' with Python...\e[0m"
conda create -y -n EGAP_env \
  "python>=3.8" \
  egap \
  -c bioconda -c conda-forge

check_success "Creating EGAP_env with Python Libraries"

# Activate the environment
echo -e "\e[36m\nActivating 'EGAP_env' environment...\e[0m"
conda activate EGAP_env
check_success "Activating 'EGAP_env' environment"

# Remove the existing runner folder if it exists
if [ -d "runner" ]; then
    echo -e "\e[33mRunner folder exists. Removing...\e[0m"
    rm -rf runner
    check_success "Removing existing runner folder"
fi

# Clone runner and install inside EGAP_env
echo -e "\e[36m\nCloning and installing runner in EGAP_env...\e[0m"
git clone https://github.com/dfguan/runner.git
check_success "Cloning runner repository"

cd runner || { echo "Failed to enter runner directory."; exit 1; }
pip install . --no-cache-dir
check_success "Installing runner via pip"

cd ..
rm -rf runner
check_success "Cleaning up runner repository"

# Verify installation
python -c "import runner" 2>/dev/null && echo -e "\e[32mRunner installed successfully in EGAP_env.\e[0m" || echo -e "\e[31mRunner installation failed.\e[0m"

# Execute quast-download commands
echo -e "\e[36m\nDownloading full suite of quast tools...\e[0m"
quast-download-gridss && quast-download-silva
check_success "Downloading quast tools"

echo -e "\n\e[32mEGAP Pipeline Pre-requisites have been successfully installed!\e[0m\n"

echo "The command 'EGAP' has been installed in your EGAP_env environment."
echo -e '\n\e[32mStart by activating the environment with "conda activate EGAP_env"\e[0m\n'
