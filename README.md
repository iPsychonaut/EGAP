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

## Setup
It is recommended that you create a dedicated conda environment for the pipeline.
```bash
conda create -y --name egap_env python=3.9 && conda activate egap_env
```

Run the main setup module.
```bash
python EGAP_setup.py
```

Command Line Example:
```bash
python EGAP_main.py -i /path/to/folder -k STRING -g INTEGER -p STRING -d [same/different] -r INTEGER -a [0/1]
```

### Command-Line Arguments
- `-i`, `--input_dir`: Path to the folder containing input data. This directory must have sub-folders named either 'illumina' (containing raw PE150 .fq.gz files and their matching MD5.txt file) or 'ont' (containing raw pass .fastq.gz files).
- `-k`, `--organism_kingdom`: Specifies the organism kingdom. Valid values are: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
- `-g`, `--genome_size`: Expected genome size in Mega-Bytes/Bases.
- `-`, `--primer_type`: String representing the Illumina primer type to use with trimmomatic. Example: 'TruSeq3-PE'.
- `-d`, `--organism_data`: Specifies whether the organism data is the same or different.
- `-r`, `--resource_use`: Integer value specifying resource usage.
- `-a`, `--attempt_install`: Flag to indicate if the installation should be attempted 0 = No (default); 1 = Yes.

## Main Functionalities
- `get_env_dir`: Retrieves environment-related directories.
- `PILON_POLISH_PIPELINE`: Polishing pipeline using PILON.
- `SPADES_HYBRID_PIPELINE`: Hybrid assembly pipeline using SPADES.

## Dependencies
Necessary Python LIbraries for this module include:
- mamba v1.5.0
- gdown v4.7.1
- biopython v1.81
- tqdm v4.38.0
- psutil v5.9.5
- termcolor v2.3.0
- beautifulsoup4 v4.12.2
- fastqc v0.11.8
- quast v5.2.0
- nanostat v1.6.0
- flye v2.9.2
- bbtools v37.62
- metaeuk v6.a5d39d9
- NanoStat v1.6.0
- minimap2 v2.26
- samtools v1.17
- bwa v0.7.17
- busco v5.5.0
- blast v2.14.1
  - makeblastdb
  - blastn
  - tblastn
- python v3.10
  - pandas
  - sys
  - subprocess
  - argparse
  - multiprocessing
  - math
  - os
  - platform
  - shutil
 
Future Implementation Planned
- abyss v2.3.7 (not used yet)
- spades.py v3.15.3 (not used yet)
- racon v1.5.0 (not used yet)
- nanoq v0.10.0 (not used yet)

Necessary JAR Files for this module include (* = version):
- trimmomatic-*.jar v0.39
- pilon*.jar v1.24

Running the main pipeline (EGAP_main.py) will attempt to install all of these prerequisites and jar files. This can be turned off by adjusting the default_install variable in EGAP_main.py.

```bash
default_install = '0' # 0 = Skip install; 1 = Attempt to install
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
The EGAP pipeline has several software prerequisites. While the pipeline attempts to install them as part of the main execution process, users are advised to ensure that all required software tools are available on their system. The `EGAP_setup` module aids in verifying and managing the presence of these prerequisites.

### EGAP_setup.py

#### Description
This module is responsible for setting up the initial environment for the EGAP pipeline. It performs tasks such as:
- Importing essential Python libraries like `os`, `re`, `tarfile`, `zipfile`, and `subprocess`.
- Using `ThreadPoolExecutor` for possible parallel execution.
- Ensuring that `mamba==1.5.0` is installed via a conda command.

#### Usage
To execute this module, run the following command:
```
python EGAP_setup.py
```


### EGAP Main Pipeline

#### Description
This module serves as the main entry point for executing the entire EGAP pipeline. It performs tasks such as:
- Providing a comprehensive command-line interface with various options like `--input_dir`, `-k`, `-g`, `-p`, `-d`, `-r`, and `-a`.
- Requiring the `--input_dir` to contain specific sub-folders for Illumina and ONT data.
- Orchestrating the execution of various steps and modules in the EGAP pipeline.

#### Usage
To execute this module, use the following command-line example:
```
python EGAP_main.py -i /path/to/folder -k STRING -g INTEGER -p STRING -d [same/different] -r INTEGER -a [0/1]
```


### EGAP Tool Checker Module

#### Description
This module performs various checks and setups to ensure that the EGAP pipeline can run smoothly. It is responsible for:
- Importing a wide range of Python libraries, including system and third-party libraries.
- Utilizing `ThreadPoolExecutor` for potential parallel execution.
- A function named `get_resource_values` to convert user-defined resource percentages into usable CPU threads and RAM values.

#### Usage
To execute this module, run the following command:
```
python check_tools.py
```


### EGAP Logging Module

#### Description
This utility module is focused on generating and managing log files within the EGAP pipeline. It performs tasks such as:
- Importing Python libraries like `pathlib` and `os` for file and directory management.
- Utilizing `termcolor` for colored terminal text and `datetime` for date and time operations.
- Providing a function named `generate_log_file` that either creates a new log file or clears an existing one.
- Another function `log_print` that prints messages to both the console and a log file.

#### Usage
To use the logging features in other modules, import this module and use its functions as needed.

```python
from log_print import log_print, generate_log_file

# Generating a log file
generate_log_file("/path/to/log_file.log")

# Logging a message
log_print("This is a message", "info")
```


### EGAP Cleaner Module

#### Description
This module focuses on cleaning up FASTA files and preparing them for further processing in the EGAP pipeline. It performs tasks such as:
- Importing essential Python libraries like `os`, `subprocess`, `glob`, `random`, `tempfile`, and `shutil`.
- Utilizing `pandas` for possible data manipulation.
- Providing a command-line interface for specifying `--dirty_fasta`, `--output_dir`, and `--organism_kingdom`.

#### Usage
To execute this module, use the following command-line example:
```
python EGAP_cleaner.py --dirty_fasta /path/to/dirty.fasta --output_dir /path/to/folder --organism_kingdom STRING
```
The `--organism_kingdom` must be from the following: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.


### EGAP Illumina Module

#### Description
This module is designed for handling Illumina sequencing data within the EGAP pipeline. It performs tasks such as:
- Providing a command-line interface for specifying `--illumina_folder`, `--primer_type`, and `-r` (an integer parameter).
- Requiring the `--illumina_folder` to contain raw Illumina `.fq.gz` files and their matching MD5.txt file.
- Requiring the `--primer_type` to be a string that represents the Illumina primer type to use with trimmomatic.

#### Usage
To execute this module, use the following command-line example:
```
python EGAP_illumina.py -i /path/to/folder -p STRING -r INTEGER
```


### EGAP ONT Module

#### Description
This module is designed for handling Oxford Nanopore Technologies (ONT) sequencing data within the EGAP pipeline. It performs tasks such as:
- Providing a command-line interface for specifying `--input_dir`, `-k`, `-g`, and `-r`.
- Requiring the `-k`, `--organism_kingdom` to be one of the following: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
- Importing a range of Python libraries for tasks like file handling, subprocess execution, and multithreading.

#### Usage
To execute this module, use the following command-line example:
```
python EGAP_ONT.py -i /path/to/folder -k STRING -g INTEGER -r INTEGER
```


### EGAP Pilon Polish Module

#### Description
This module is designed for genome polishing using Pilon. It typically uses both long-read (like ONT) and short-read (like Illumina) sequencing data for high-quality genome assembly. It performs tasks such as:
- Providing a command-line interface for specifying `-oi`, `-if`, `-ir`, `-k`, and `-r`.
- Requiring the `-k`, `--organism_kingdom` to be one of the following: Archaea, Bacteria, Fauna, Flora, Funga, or Protista.
- Importing a range of Python libraries for tasks like file handling, subprocess execution, and multithreading.

#### Usage
To execute this module, use the following command-line example:
```
python EGAP_pilon_polish.py -oi /path/to/cleaned_ont_reads.fastq -if /path/to/forward_reads.fq -ir /path/to/reverse_reads.fq -k STRING -r INTEGER
```


### EGAP Quality Control (QC) Module

#### Description
This module is focused on quality control (QC) and assessment of genome assemblies within the EGAP pipeline. It performs tasks such as:
- Importing a range of Python libraries like `BeautifulSoup` for parsing HTML/XML documents and `pandas` for data manipulation.
- Providing a function named `assess_with_quast` to run QUAST on an assembly for quality assessment.
- Utilizing custom functions like `log_print` and `generate_log_file` for logging, and `get_env_dir` for environment setup.

#### Usage
To execute this module, run the following command:
```
python EGAP_qc.py
```



## Future Expansions
- MAIN TODO: Integration of Kelsey's Slot Lab 'different organism' SPAdes assembly Racon polishing pipeline. This includes getting abyss-sealer working.
- FASTQC TODO: If warning is in Overrepresented sequences AND?/OR? Adapter Content then there is likely a Failed Primers Trimming Error.
- OVERALL TODO: Perform file removal to reduce data bloat.



## License
EGAP is distributed under the BSD 3-Clause License.