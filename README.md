
# EGAP Main Module

## Description
This module serves as the primary entry point for the EGAP program. It provides functionalities for genome assembly and polishing using various tools and pipelines.

## Author Information
- **Name**: Ian Michael Bollinger
- **Email**: ian.michael.bollinger@gmail.com
- **Collaboration**: Developed with the help of ChatGPT 4.0

## Intended Platform
This program is intended and tested for WSL/Linux based systems. It is recommended to use it on such platforms for optimal performance and compatibility.

## Installation

### Directory Structure
The program consists of eight (8) modules and two (2) database folders. It is recommended to save all these components in a base-level directory of the current operating system named 'EGAP'.

### Cloning the Repository
To clone the entire EGAP repository to the desired location and extract out the necessary databases for processing, use the following command:
```bash
git clone https://github.com/iPsychonaut/EGAP.git ~/EGAP
```

## Setup and Installation

### Initial Setup with Bash Scripts

To facilitate an easy and robust setup process, we have divided the setup into two parts, each handled by a separate Bash script. This approach ensures a smoother installation and setup experience, especially in managing dependencies and environment configurations.

#### Part 1: Preparing the Environment

1. Run the `EGAP_setup.sh` script. This script handles the initial setup, including:
   - Updating the Java Runtime Environment (JRE).
   - Installing unzip.
   - Checking and updating Conda.
   - Creating a custom Conda environment named 'entheome_env'.

   To run the script, use the following command in your terminal:
   ```bash
   bash EGAP_setup.sh
   ```

2. After successfully running `EGAP_setup.sh`, **close and restart your terminal** to ensure that all changes are properly applied.

#### Part 2: Completing the Setup

1. Once you have restarted your terminal, activate the entheome_env.
   ```bash
   conda activate entheome_env
   ```

2. Now install all necessary Python libraries with Mamba:
   ```bash
   mamba install -y -c bioconda -c conda-forge -c agbiome -c prkrekel pandas==2.0.3 busco==5.5.0 nanoq==0.10.0 biopython==1.81 tqdm==4.38.0 beautifulsoup4==4.12.2 quast==5.2.0 nanostat==1.6.0 flye==2.9.2 bbtools==37.62 metaeuk==6.a5d39d9 blast==2.14.1 bwa==0.7.17 minimap2==2.26 pysam==0.21.0 samtools==1.17 arcs==1.2.5 tigmint==1.2.10 abyss==2.3.7 racon==1.5.0 spades==3.15.3 gdown==4.7.1 psutil==5.9.5 abyss==2.3.7 requests==2.31.0 minimap2==2.26 spoa==4.1.3 racon==1.5.0 termcolor==2.3.0 fastqc==0.12.1 masurca==4.1.0 openjdk=8
   ```

### Note on Script Execution

- Ensure that you have the necessary permissions to execute these scripts. You might need to run `chmod +x EGAP_setup_1.sh EGAP_setup_2.sh` to make the scripts executable.
- It's crucial to follow the order of the scripts and restart your terminal between executing `EGAP_setup_1.sh` and `EGAP_setup_2.sh` for proper setup.
- These scripts are designed to check for errors at each step and will stop execution if an error is encountered. This is to prevent cascading errors and ensure a reliable setup process.

### Post-Setup Actions

After successfully running both scripts, your environment should be ready for the EGAP pipeline. To activate the environment, use:
```bash
conda activate entheome_env
```

Now, you are all set to proceed with using the EGAP program for genome assembly and analysis!
Command Line Example:
```bash
python EGAP_main.py --input_dir /path/to/folder --reads_type hybrid --organism_kingdom STRING --genome_size INTEGER --primer_type STRING --org_data same/different --resource_use INTEGER
```

### Command-Line Arguments
- `--input_dir`: Path to the folder containing input data. This directory must have sub-folders named either 'illumina' (containing raw PE150 .fq.gz files and their matching MD5.txt file) or 'ont' (containing raw pass .fastq.gz files).
- `--reads_data`: 
- `--organism_kingdom`: Specifies the organism kingdom. Valid values are: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
- `--genome_size`: Expected genome size in Mega-Bytes/Bases.
- `--primer_type`: String representing the Illumina primer type to use with trimmomatic. Example: 'TruSeq3-PE'.
- `--org_data`: Specifies whether the organism data is the same or different.
- `--resource_use`: Integer value specifying resource usage.

## Main Functionalities
- `get_env_dir`: Retrieves environment-related directories.
- `PILON_POLISH_PIPELINE`: Hybrid Assembly Polishing pipeline using PILON.

## Dependencies
Necessary Python LIbraries for this module include:
- pandas
- biopython
- tqdm
- psutil
- beautifulsoup4
- fastqc
- quast
- nanoq
- nanostat
- flye
- bbtools
- metaeuk
- sys
- subprocess
- argparse
- multiprocessing
- math
- os
- platform
- shutil

Necessary Third-Party programs for this module include:
- busco
- compleasm
- samtools
- bwa
- makeblastdb
- blastn
- racon (not used yet)
- nanoq
- NanoStat
- spades.py (not used yet)
- tblastn
- flye
- minimap2
- metaeuk
- tqdm

Necessary JAR Files for this module include (* = version):
- trimmomatic-*.jar
- pilon*.jar

Running the main pipeline (EGAP_main.py) will attempt to install all of these prerequisites and jar files. This can be turned off by adjusting the default_install variable in EGAP_main.py.

```bash
default_install = '0' # 0 = Attempt to install; 1 = Skip install
```

## Contact Information
For any queries or issues related to the EGAP program, please contact:
- **Email**: ian.michael.bollinger@gmail.com


## Overview of the EGAP Pipeline

EGAP (Example Genome Assembly Pipeline) is a comprehensive tool designed to process and analyze genomic data. The pipeline is modular, consisting of several Python modules, each dedicated to a specific task in the processing chain.

### Standalone Modules
The following modules are designed to run independently:
- `EGAP_illumina`: Processes Illumina reads.
- `EGAP_ONT`: Processes Oxford Nanopore Technology (ONT) reads.
- `EGAP_pilon_polish`: Polishes assemblies using Pilon.

### Integrated Pipeline
For a streamlined experience, users can execute the `EGAP_main` module, which will orchestrate the entire analysis pipeline, integrating functionalities from all modules.

### Prerequisites
The EGAP pipeline has several software prerequisites. While the pipeline attempts to install them as part of the main execution process, users are advised to ensure that all required software tools are available on their system. The `check_tools` module aids in verifying and managing the presence of these prerequisites.


## EGAP Illumina Module

### Description
This module handles various operations related to Illumina data processing, including extraction, verification, trimming, and overall processing. Currently the pipeline can determine if the data is Single-End (SE) or Paired-END (PE); however it has only been tested with PE data; if you have SE data that can be tested please contact the author.

### Main Functionalities
- `md5_check`: Verifies the integrity of files using MD5 hashing.
- `illumina_extract_and_check`: Manages the extraction and verification of Illumina data.
- `trim_with_trimmomatic`: Trims Illumina reads using the Trimmomatic tool.
- `process_illumina`: Orchestrates the processing of Illumina data.


## EGAP ONT Module

### Description
This module focuses on the processing and assembly of Oxford Nanopore Technologies (ONT) sequencing data. It provides tools for combining raw data, assembling sequences using Flye, and overall ONT data management.

### Main Functionalities
- `ont_combine_fastq_gz`: Combines multiple fastq.gz files related to ONT data.
- `assemble_ont_flye`: Assembles sequences using the Flye tool tailored for ONT data.
- `process_ONT`: Manages the overall processing of ONT data.


## EGAP Cleaner Module

### Description
This module is dedicated to cleaning and managing sequence data. It provides tools for creating BLAST databases, handling sequence chunks, running BLASTn operations, and cleaning FASTA files. Provided are model organisms (Assembled_Databases) for each of the major kingdoms, however this is not an extensive list; feel free to include other organisms of interest.

### Main Functionalities
- `create_BLAST_db`: Creates a BLAST database for sequences.
- `chunk_seq`: Divides sequences into manageable chunks or segments.
- `handle_chunk`: Manages individual sequence chunks.
- `run_blastn`: Conducts BLASTn operations on sequences.
- `clean_dirty_fasta`: Cleans or sanitizes FASTA files to remove unwanted sequences or contaminants.


## EGAP Pilon Polish Module

### Description
This module focuses on the polishing of ONT Flye assemblies using Illumina reads. It integrates tools like BWA, Samtools, and Pilon to achieve a polished assembly from raw reads.

### Main Functionalities
- `index_cleaned_ont_assembly`: Indexes a cleaned ONT Flye assembly using BWA.
- `generate_illumina_sam`: Generates a SAM map of the cleaned ONT Flye de novo assembly using the trimmed Illumina paired forward & reverse reads.
- `convert_sam_to_bam`: Converts a SAM file into a BAM file using Samtools.
- `pilon_polish_assembly`: Polishes a cleaned ONT Flye assembly using Illumina Binary Alignment Map with Pilon.
- `final_pilon_polish`: Main function that handles the polishing of cleaned ONT reads with trimmed Illumina reads.


## EGAP Quality Control (QC) Module

### Description
This module is dedicated to assessing the quality of sequencing data and assemblies. It integrates various quality control tools like QUAST, BUSCO, FastQC, and NanoStat.

### Main Functionalities
- `assess_with_quast`: Conducts an assessment of assemblies using QUAST.
- `assess_with_busco`: Evaluates assemblies for the presence of universal single-copy orthologs using BUSCO.
- `assess_with_fastqc`: Runs a quality assessment on raw sequencing reads using FastQC.
- `assess_with_nanostat`: Provides statistics for Oxford Nanopore sequencing data using NanoStat.


## EGAP Tool Checker Module

### Description
This module provides functionalities for checking, finding, and managing files and software prerequisites necessary for the operation of the main program.

### Main Functionalities
- `move_file_up`: Moves a specified file up in the directory hierarchy.
- `get_md5`: Computes the MD5 hash for a given file.
- `search_directory_for_file`: Searches within a directory for a specific file.
- `find_file`: Finds a specific file within a directory or its subdirectories.
- `get_env_dir`: Retrieves the directory of the current environment.
- `libraries_check`: Checks for the presence of specific libraries.
- `check_for_jars`: Checks for the presence of specific Java JAR files.
- `check_prereqs_installed`: Checks if prerequisite software is installed.


## EGAP Logging Module

### Description
This module is dedicated to logging messages. It aids in recording progress, results, and any messages during the execution of the program. It offers dual logging capabilities - to the console and to a specified log file.

### Main Functionalities
- `generate_log_file`: Generates a log file to record progress and results.
- `log_print`: Logs messages simultaneously to the console and the specified log file.


## Future Expansions
- MAIN TODO: Add versions for each Prerequisite/Third-Party program as well as python version that maintains functionality.
- MAIN TODO: Integration of Kelsey's Slot Lab 'different organism' SPAdes assembly Racon polishing pipeline. This includes getting abyss-sealer working.
- FASTQC TODO: If warning is in Overrepresented sequences AND?/OR? Adapter Content then there is likely a Failed Primers Trimming Error.

## License
EGAP is distributed under the BSD 3-Clause License.
