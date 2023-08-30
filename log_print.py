# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 08:53:45 2023

@author: ian.michael.bollinger@gmail.com with the help of ChatGPT 4.0

"""
import pathlib, os
from termcolor import cprint
from datetime import datetime

# Function to Generate a Log File
def generate_log_file(log_file_path, use_numerical_suffix=False):
    """
    Generates or clears a log file based on the given parameters.
    
    This function either creates a new log file or clears an existing one, depending on the specified parameters. 
    If the 'use_numerical_suffix' parameter is True, and the file already exists, a new file with a numerical suffix 
    will be created. Otherwise, the existing file will be cleared.
    
    Parameters:
        log_file_path (str): Path to the log file.
        use_numerical_suffix (bool): If True, creates new files with numerical suffixes if the file exists; otherwise, clears the existing file.
    
    Returns:
        str: Path to the log file.
    """
    if os.path.exists(log_file_path) and use_numerical_suffix:
        # If using numerical suffixes, increment until a new filename is found
        counter = 1
        new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        while os.path.exists(new_log_file_path):
            counter += 1
            new_log_file_path = f"{log_file_path.rsplit('.', 1)[0]}_{counter}.txt"
        log_file_path = new_log_file_path
    else:
        # Clear the existing log file or create a new one
        open(log_file_path, 'w').close()
    
    return log_file_path

# Function to color coded print to console and save to log_file information
def log_print(input_message, log_file):
    """
    Logs a message to a file and prints it to the console with appropriate coloring.
    
    This function takes a message and logs it to the specified file. Additionally, the message is printed to the 
    console, potentially with specific coloring depending on the context.
    
    Parameters:
        input_message (str): Message to be logged and printed.
        log_file (str): Path to the log file.
    """
    # Establish current date-time
    now = datetime.now()
    message = f'[{now:%Y-%m-%d %H:%M:%S}]\t{input_message}'

    # Determine the print color based on the input_message content
    message_type_dict = {'ERROR': ['red'],
                         'WARN': ['yellow'],
                         'NOTE': ['blue'],
                         'CMD': ['cyan'],
                         'PASS': ['green']}
    print_color = ['white']  # Default color
    for key, value in message_type_dict.items():
        if key.lower() in input_message.lower():
            print_color = value
            break

    # Writing the message to the log file
    with open(log_file, 'a') as file:
        print(message, file=file)

    # Handling different message types for colored printing
    try:
        cprint(message, print_color[0])
    except (KeyError, IndexError):
        cprint(message, print_color[1] if len(print_color) > 1 else 'white')

## Debuging Main Space & Example
if __name__ == '__main__':
    print('log_print.py DEBUG')
    
    # Establish Environment
    input_file = 'E:/EGAP/test.fastq'
    try:
        with open(input_file, 'rb') as f:
            print('WINDOWS ENVIRONMENT')
    except FileNotFoundError:
        print('LINUX/WSL/MAC ENVIRONMENT')
        input_file = input_file.replace(input_file.split('/')[0],'/mnt/e')    

    # Generate log file with the desired behavior
    file_extension = pathlib.Path(input_file).suffix
    run_log = input_file.replace(file_extension, '.tsv')
    run_log = generate_log_file(run_log, use_numerical_suffix=False)

    log_print(input_message='PASS: TEST', log_file=run_log)
    log_print(input_message='NOTE: TEST', log_file=run_log)
    log_print(input_message='ERROR: TEST', log_file=run_log)
    log_print(input_message='WARN: TEST', log_file=run_log)
    log_print(input_message='CMD: TEST', log_file=run_log)
    log_print(input_message='TEST', log_file=run_log)
