
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

1. Run the `EGAP_setup_1.sh` script. This script handles the initial setup, including:
   - Updating the Java Runtime Environment (JRE).
   - Installing unzip.
   - Checking and updating Conda.
   - Creating a custom Conda environment named 'entheome_env'.

   To run the script, use the following command in your terminal:
   ```bash
   bash EGAP_setup_1.sh
   ```

2. After successfully running `EGAP_setup_1.sh`, **close and restart your terminal** to ensure that all changes are properly applied.

#### Part 2: Completing the Setup

1. Once you have restarted your terminal, run the `EGAP_setup_2.sh` script. This script completes the setup process by:
   - Continuing with the installation of necessary tools and libraries.
   - Setting up additional environment configurations.

   Execute the script with the command:
   ```bash
   bash EGAP_setup_2.sh
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

[... Additional content of the original README ...]

## License
EGAP is distributed under the BSD 3-Clause License.
