#!/bin/bash

# Perform JRE update and check
echo -e "\e[36mUpdating Java Runtime Environment (JRE)...\e[0m"
sudo apt-get update
sudo apt install openjdk-8-jre-headless -y

# While in Sudo Get Unzip
echo -e "\e[36mDownloading unzip...\e[0m"
sudo apt install unzip

# Update pip in the base environment
echo -e "\e[36mUpdating pip...\e[0m"
conda run -n base pip install --upgrade pip

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo -e "\e[36mConda is not installed. Please install Conda before running this script."
    echo "Run the following commands to install Miniconda3:"
    echo "     wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "     bash Miniconda3-latest-Linux-x86_64.sh"
    echo "     source ~/.bashrc"
    exit 1
else
    # Update Conda
    echo -e "\e[36mUpdating Conda...\e[0m"
    conda update -n base -c defaults conda -y
fi

# Install and update Mamba in the base environment
echo -e "\e[36mInstalling and updating Mamba...\e[0m"
conda install mamba -n base -c conda-forge -y
conda update mamba -n base -c conda-forge -y

# Create a custom Conda environment with Python 3.8.15 for main Entheome Pipeline
echo -e "\e[36mCreating custom Conda environment 'entheome_env' with Python 3.8.15...\e[0m"
conda create -n entheome_env python=3.8.15 -y

# Check if gdown is installed, and install it if it's not
if ! command -v gdown &> /dev/null; then
    echo -e "\e[36mInstalling gdown...\e[0m"
    conda run -n entheome_env pip install gdown==4.7.1
fi

# Download the Google Drive file
echo -e "\e[36mDownloading databases from Google Drive...\e[0m"
conda run -n entheome_env gdown https://drive.google.com/uc?id=1Hj-8tFlJPiOoP_8_Sp4pyWipaQ3zr745 -O EGAP_Databases.zip

# Create the EGAP_Databases directory and extract the contents of the zip file
echo -e "\e[36mExtracting contents to './EGAP_Databases'...\e[0m"
mkdir -p ./EGAP_Databases
unzip EGAP_Databases.zip -d ./EGAP_Databases
rm EGAP_Databases.zip

# Make Sure all Mamba installations are complete
echo -e "\e[36mDownloading Python libraries...\e[0m"
conda run -n entheome_env mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel pandas==2.0.3 busco==5.5.0 nanoq==0.10.0 biopython==1.81 tqdm==4.38.0 beautifulsoup4==4.12.2 fastqc==0.11.8 quast==5.2.0 nanostat==1.6.0 flye==2.9.2 bbtools==37.62 metaeuk==6.a5d39d9 blast==2.14.1 bwa==0.7.17 minimap2==2.26 pysam==0.21.0 samtools==1.17 arcs==1.2.5 tigmint==1.2.10 abyss==2.3.7 racon==1.5.0 spades==3.15.3 gdown==4.7.1 psutil==5.9.5 abyss==2.3.7 requests==2.31.0 minimap2==2.26 spoa==4.1.3 racon==1.5.0 termcolor==2.3.0 fastqc==0.12.1 masurca==4.1.0 openjdk=8

# sepp==4.5.1 into funannotate env

# Download and set up Trimmomatic
echo -e "\e[36mDownloading and setting up Trimmomatic...\e[0m"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm ./Trimmomatic-0.39.zip
# make ./Trimmomatic-0.39/Trimmomatic-0.39.jar callable as Trimmomatic-0.39.jar


# Download and set up Pilon
echo -e "\e[36mDownloading and setting up Pilon...\e[0m"
wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar
mkdir -p Pilon
mv pilon-1.24.jar Pilon/

# Download and set up Compleasm
echo -e "\e[36mDownloading and setting up Compleasm...\e[0m"
wget https://github.com/huangnengCSU/compleasm/releases/download/v0.2.2/compleasm-0.2.2_x64-linux.tar.bz2
mkdir -p Compleasm
tar -jxvf compleasm-0.2.2_x64-linux.tar.bz2 -C Compleasm
chmod +x Compleasm/compleasm_kit/compleasm.py
rm ./compleasm-0.2.2_x64-linux.tar.bz2

# Get the absolute path of the directory where this script is located
# echo -e "\e[36mUpdating ~/.bashrc...\e[0m"
# SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# Add exports to .bashrc for future use
# echo "export GENEMARK_PATH=$SCRIPT_DIR/GeneMark/gmes_linux_64_4" >> ~/.bashrc
# echo "export FUNANNOTATE_DB=$SCRIPT_DIR/EGAP_Databases/Funannotate_Databases" >> ~/.bashrc

# Add aliases to .bashrc for easy access
# echo "alias pilon='java -jar $SCRIPT_DIR/Pilon/pilon-1.24.jar'" >> ~/.bashrc
# echo "alias compleasm='python $SCRIPT_DIR/Compleasm/compleasm_kit/compleasm.py'" >> ~/.bashrc

# Source .bashrc to apply the changes
echo -e "\e[33mPlease run 'conda activate entheome_env' to begin processing.\e[0m"
echo -e "\e[32mEntheome pipeline set-up complete.\e[0m"
