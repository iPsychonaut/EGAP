# -*- coding: utf-8 -*-
"""
Created on Sun Aug 20 12:22:32 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0
"""
import os, subprocess

# Ensure all other libraries are installed
libraries = ['pandas', 'biopython', 'tqdm', 'psutil',
             'termcolor', 'beautifulsoup4', 'fastqc',
             "quast", 'nanoq', 'nanostat', 'flye']
print(f'UNLOGGED:\tAttempting to install the following Python Libraries: {libraries}')
try:
    # Run a subprocess calling conda install for each missing library
    exit_code = subprocess.check_call(['conda', 'install', '-y', '-c', 'bioconda', *libraries],
                                      stdout=subprocess.DEVNULL, 
                                      stderr=subprocess.DEVNULL)
    if exit_code == 0:
        # Additional commands for downloading databases/tools for QUAST
        additional_commands = ['quast-download-gridss', 'quast-download-silva', 'quast-download-busco']
        for cmd in additional_commands:
            try:
                subprocess.check_call(cmd, shell=True,
                                      stdout=subprocess.DEVNULL, 
                                      stderr=subprocess.DEVNULL)
                print(f'UNLOGGED PASS:\tSuccessfully executed {cmd}')
            except subprocess.CalledProcessError:
                print(f"UNLOGGED ERROR:\tFailed to execute {cmd}")
    else:
        print(f"UNLOGGED ERROR:\t Unable to Install {libraries}")
except:
    # Print an error if something goes wrong during the installation
    print(f"UNLOGGED ERROR:\t Unable to Install {libraries} library")