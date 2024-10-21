# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:27:53 2023

@author: ian.michael.bollinger@gmail.com/researchconsultants@critical.consulting

Command Line Example:
    python EGAP.py -i /path/to/base/folder -d READS_DATA_STRING
                   -k ORGANISM_KINGDOM_STRING  -es ESTIMATED_GENOME_SIZE_MBP
                   -r PERCENT_RESOURCES_FLOAT  -rf /path/to/reference_sequence

Arguments:
    -i, --input_folder: Path to the input folder containing FASTQ files.
                        (default: /mnt/d/ENTHEOME/Ps_semilanceata/Illumina_PE150)
    -d, --reads_data: Specify the type of sequencing data. Must be one of:
                      'illu' for Illumina, 'ont' for Nanopore, or 'hybrid' for combined data.
                      (default: hybrid)
    -k, --organism_kingdom: The kingdom to which the current organism belongs.
                            Options: Archaea, Bacteria, Fauna, Flora, Funga, Protista.
                            (default: Funga)
    -es, --est_size: Estimated genome size in Megabase pairs (Mbp). For example, '60' will
                     be converted to 60,000,000 base pairs. (default: 60)
    -r, --resources: Percentage of available system resources (CPU and memory) to use.
                     Must be between 0.01 and 1.00. (default: 0.45)
    -rf, --ref_seq: Optional path to the reference genome for guiding the assembly.
                    (default: None)
"""
# Base Python Imports
import os, argparse


# Custom Python Imports
from EGAP_Tools import log_print, initialize_logging_environment
from EGAP_Illumina import illumina_prep, illumina_only_main
from EGAP_ONT import ont_prep, ont_only_main
from EGAP_Masurca import masurca_de_novo, masurca_ref_seq
from EGAP_Pilon import final_hybrid_pilon
from EGAP_QC import assembly_qc_checks


# Global output_area variable
CPU_THREADS = 1
RAM_GB = 1
DEFAULT_LOG_FILE = None
ENVIRONMENT_TYPE = None


## Debuging Main Space & Example
if __name__ == "__main__":    
    # Define and parse command line arguments
    parser = argparse.ArgumentParser(description="Run Entheome Genome Assembly Pipeline")

    # Default values
    default_input_folder = '/mnt/d/ENTHEOME/Ps_semilanceata/'
    default_reads_data = 'hybrid'
    default_organism_kingdom = 'Funga'
    default_estimated_genome_size = 60
    default_percent_resources = 0.75
    default_reference_sequence = None

    # Add arguments
    parser.add_argument("-i", "--input_folder", type=str, default=default_input_folder,
                        help="Path to the input folder containing FASTQ files")
    parser.add_argument('--reads_data', '-d',
                        type = str, default = default_reads_data,
                        choices = ["illu","ont","hybrid"],
                        help = f'Indicate if the provided data are generated from the same organism or different organisms (default: {default_reads_data})')
    parser.add_argument('--organism_kingdom', '-k',
                        type=str, default=default_organism_kingdom,
                        help = f'Kingdom the current organism data belongs to. (default: {default_organism_kingdom})')
    parser.add_argument("--est_size", "-es",
                        type=int, default=default_estimated_genome_size,
                        help="Estimaged size of the genome in Mbp (aka million-base-pairs). (default: {default_estimated_genome_size}Mbp => {int(default_estimated_genome_size)*1000000})")
    parser.add_argument("-r", "--resources",
                        type=float, default=default_percent_resources,
                        help=f"Percentage of resources to use. (0.01-1.00; default: {default_percent_resources})")
    parser.add_argument("--ref_seq", "-rf",
                        type=str, default=default_reference_sequence,
                        help="Path to the reference genome for assembly. (default: {default_reference_sequence})")

    # Parse the arguments
    args = parser.parse_args()
    BASE_FOLDER = args.input_folder
    READS_DATA = args.reads_data
    CURRENT_ORGANISM_KINGDOM = args.organism_kingdom
    EST_SIZE = args.est_size * 1000000
    PERCENT_RESOURCES = args.resources
    REF_SEQ = args.ref_seq
    
    # Generate log file with the desired behavior
    initialize_logging_environment(BASE_FOLDER)
        
    # Run appropriate pipe according READS_DATA
    if READS_DATA == 'illu':
        log_print("RUN ILLUMINA ONLY PIPELINE")

        masurca_scaffolded_assmebly_path = illumina_only_main(BASE_FOLDER,
                                                              READS_DATA,
                                                              CURRENT_ORGANISM_KINGDOM,
                                                              EST_SIZE,
                                                              PERCENT_RESOURCES,
                                                              REF_SEQ)

        log_print("PASS:\Illumina Assembly Sucessfully Polished; EGAP Illumina Assembly Complete!")

    elif READS_DATA == 'ont':
        log_print("RUN ONT ONLY PIPELINE")

        output_assembly_path = ont_only_main(BASE_FOLDER, READS_DATA,
                                             CURRENT_ORGANISM_KINGDOM,
                                             EST_SIZE, PERCENT_RESOURCES, REF_SEQ)

        log_print("PASS:\ONT Assembly Sucessfully Polished; EGAP ONT Assembly Complete!")

    elif READS_DATA == 'hybrid':
        log_print("RUN HYBRID PIPELINE")
        
        ILLU_FOLDER = os.path.join(BASE_FOLDER, "Illumina_PE150")

        ONT_FOLDER = os.path.join(BASE_FOLDER, "ONT_MinION")
        
        input_fq_list, fastqc_output_dirs = illumina_prep(ILLU_FOLDER, READS_DATA,
                                                          CURRENT_ORGANISM_KINGDOM,
                                                          EST_SIZE, PERCENT_RESOURCES,
                                                          CPU_THREADS, RAM_GB, REF_SEQ)
        
        print(input_fq_list)
        
        combined_ont_fastq = ont_prep(ONT_FOLDER, READS_DATA,
                                      CURRENT_ORGANISM_KINGDOM, EST_SIZE,
                                      PERCENT_RESOURCES, REF_SEQ)
        de_novo_output_folder = os.path.join(BASE_FOLDER, "Hybrid-De-Novo-Assembly")
        ref_seq_output_folder = os.path.join(BASE_FOLDER, "Ref_Seq_Assembly")

        if REF_SEQ != None:        

            final_ref_seq_assembly_path, bam_list = masurca_ref_seq(BASE_FOLDER, de_novo_output_folder,
                                                                    [input_fq_list[0], input_fq_list[1], combined_ont_fastq],
                                                                    EST_SIZE, PERCENT_RESOURCES, REF_SEQ)
    
            log_print("PASS:\tIllumina & ONT reads Sucessfully Reference Sequence Assembled; ready for next steps.\n")

            pilon_polished_assembly = final_hybrid_pilon(final_ref_seq_assembly_path,
                                                         [input_fq_list[0], input_fq_list[1]],
                                                         READS_DATA, CURRENT_ORGANISM_KINGDOM,
                                                         PERCENT_RESOURCES, REF_SEQ)
        else:
            
            final_de_novo_assembly_path, bam_list = masurca_de_novo(BASE_FOLDER, de_novo_output_folder,
                                                                    [input_fq_list[0], input_fq_list[1], combined_ont_fastq],
                                                                    EST_SIZE, PERCENT_RESOURCES, REF_SEQ)
    
            log_print("PASS:\tIllumina & ONT reads Sucessfully de novo Assembled; ready for next steps.\n")

            pilon_polished_assembly = final_hybrid_pilon(final_de_novo_assembly_path,
                                                         [input_fq_list[0], input_fq_list[1]],
                                                         READS_DATA, CURRENT_ORGANISM_KINGDOM,
                                                         PERCENT_RESOURCES, REF_SEQ)

        # Quality Control Check Assembly with QUAST & Compleasm
        quast_output_dir, compleasm_output_dir_1, compleasm_output_dir_2 = assembly_qc_checks(pilon_polished_assembly, READS_DATA, CURRENT_ORGANISM_KINGDOM, REF_SEQ, CPU_THREADS)
        
        log_print("PASS:\tHybrid Assembly Sucessfully Polished; EGAP Hybrid Assembly Complete!")