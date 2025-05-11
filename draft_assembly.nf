#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//
// EGAP Version
//
def version = "3.0.0f"


//
// Define Parameters
//
params.no_file = "$projectDir/assets/NO_FILE"
params.input_csv = "/mnt/d/EGAP_Nextflow/EGAP_test.csv"
params.output_dir = "/mnt/d/TESTING_SPACE/nextflow_test"
params.cpu_threads = 12
params.ram_gb = 40


//
// Print Banner & Parameters
//
log.info("""
\033[91m.---.\033[92m________\\\033[38;5;208m=\033[96m/\033[92m________\033[91m.---.\033[0m 
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m--------\033[96m/\033[91m=\033[92m\\\033[94m--------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|   \033[94m.\033[96m---------.  \033[94m.\033[96m------.    \033[94m.\033[96m------.    \033[94m.\033[96m-------.\033[0m
\033[91m`---'\033[96m~~~~~~~(\033[38;5;208m===\033[92m)\033[96m~~~~~~~\033[91m`---'  \033[94m/\033[96m|         |\033[94m/\033[96m/        \\ \033[94m/\033[96m/        \\ \033[94m/\033[96m/         \\\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m.--     \033[96m\\\033[91m=\033[92m/    \033[92m,--. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .------\033[94m'\033[96m|   .------\033[94m'\033[96m|   .--\033[94m.   \033[96m. |   .--\033[94m.\033[96m   .\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m|-      \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|  _ \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m----"| \033[96m|   |\033[94m----"| \033[96m|   |\033[94m-'\033[96m|   | |   |\033[94m-'\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[96m`--    \033[92m(\033[91m===\033[96m)   \033[96m`--' \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `----.\033[94m|\033[96m |   |     \033[94m| \033[96m|   +--+   | |   +--+   |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[92m\\\033[38;5;208m=\033[96m/         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|       |\033[94m| \033[96m|   | \033[94m.\033[96m----.|          | |          |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[96m/\033[91m=\033[92m\\         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .----'\033[94m| \033[96m|   |\033[94m"\033[96m|    ||   +--+   | |   +------'\033[0m
 \033[92m|\033[94m|\033[96m|        \033[96m(\033[38;5;208m===\033[92m)        \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m---" | \033[96m|   |\033[94m"\033[96m`-.  ||   |\033[94m| \033[96m|   | |   |\033[94m-----"\033[0m
 \033[92m|\033[94m|\033[96m|  \033[96m__     \033[96m\\\033[91m=\033[92m/    \033[96m,__. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `------.|   `---'  ||   |\033[94m| \033[96m|   | |   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m/__\\    \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|__| \033[96m|\033[94m|\033[92m|  \033[94m`.\033[96m|         |`.        \033[96m/\033[94m.\033[96m|   |\033[94m| \033[96m|   |\033[94m.\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m|  |   \033[92m(\033[91m===\033[96m)   \033[92m|    \033[96m|\033[94m|\033[92m|    \033[96m`---------'  `------'  `---' \033[94m`\033[96m'---' `---'\033[0m
\033[91m.---.\033[96m_____\033[92m[:\033[94m:\033[96m::::\033[94m:\033[92m]\033[96m_____\033[91m.---.\033[92m    \033[92m╔═══════════════════════════════════════════╗\033[0m
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m------\033[92m|:\033[94m:\033[96m::\033[94m:\033[92m|\033[94m------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|    \033[92m║     \033[94mEnthe\033[96mome Ge\033[97mnome Assemb\033[96mly Pip\033[94meline     \033[92m║\033[0m
\033[91m`---'\033[92m~~~~~~\033[92m|::\033[94m:\033[96m:\033[94m:\033[92m|~~~~~~\033[91m`---'    \033[92m╚═══════════════════════════════════════════╝\033[0m

                    Curated & Maintained by Ian M Bollinger
                         (\033[94mian.bollinger@entheome.org)\033[0m

                              \033[92mdraft_assembly.nf\033[0m
                                version ${version}

 Input-Setup \033[94m-\033[92m>\033[0m Preprocess \033[94m-\033[92m>\033[0m Assemble \033[94m-\033[92m>\033[0m Compare \033[94m-\033[92m>\033[0m Polish \033[94m-\033[92m>\033[0m Curate \033[94m-\033[92m>\033[0m Assess""")
sleep(1500)
log.info("""
\033[91m================================================================================\033[0m

    \033[92minput-setup settings
        \033[94mstarted at                       \033[0m: ${workflow.start}
        \033[94mconfig files                     \033[0m: ${workflow.configFiles}
        \033[94mcontainer                        \033[0m: ${workflow.containerEngine}:${workflow.container}
        \033[94mRAM GB                           \033[0m: ${params.ram_gb}
        \033[94mCPU threads                      \033[0m: ${params.cpu_threads}
        \033[94minput csv                        \033[0m: ${params.input_csv}
        \033[94moutput to                        \033[0m: ${params.output_dir}

    \033[96m-----------------------------------------------------------------------\033[0m""")
sleep(250)
log.info("""
    \033[92mtrimmomatic settings
        \033[94mmode                             \033[0m: -PE
        \033[94mphred version                    \033[0m: -phred33
        \033[94milluminaclip adapter             \033[0m: /opt/conda/envs/EGAP_env/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa
        \033[94mfastaWithAdaptersEtc             \033[0m: 2
        \033[94mseed mismatches                  \033[0m: 30
        \033[94mpalindrome clip threshold        \033[0m: 10
        \033[94msimple clip threshold            \033[0m: 11
        \033[94mHEADCROP                         \033[0m: 10
        \033[94mCROP                             \033[0m: 145
        \033[94mSLIDINGWINDOW                    \033[0m: 50:25
        \033[94mMINLEN                           \033[0m: 125
        
    \033[92mbbduk settings
        \033[94mktrim                            \033[0m: -r
        \033[94mk                                \033[0m: 23
        \033[94mmink                             \033[0m: 11
        \033[94mhdist                            \033[0m: 1
        \033[94mtrimpairsevenly                  \033[0m: -tpe
        \033[94mtrimbyoverlap                    \033[0m: -tbo
        \033[94mqtrim                            \033[0m: -rl
        \033[94mtrimq                            \033[0m: 20

    \033[92mclumpify settings
        \033[94mremove duplicate reads           \033[0m: -dedupe

    \033[96m-----------------------------------------------------------------------\033[0m""")
sleep(250)
log.info("""
    \033[92mfiltlong settings
        \033[94mmin_length                       \033[0m: 1000
        \033[94mmin_mean_q                       \033[0m: 8
        \033[94mkeep_percent                     \033[0m: 90 
        \033[94mcoverage                         \033[0m: 75
        \033[94mtarget_bases                     \033[0m: estimated size (bp) * 75 (coverage)
        
    \033[92mratatosk settings
        \033[94mverbose output                   \033[0m: -v

    \033[96m-----------------------------------------------------------------------\033[0m""")
sleep(250)
log.info("""
    \033[92mmasurca settings
        \033[94mconfiguration file sections      \033[0m: DATA, PARAMETERS
        \033[94mgraph kmer size                  \033[0m: auto
        \033[94muse linking mates                \033[0m: (0 if Hybrid assembly, else 1)
        \033[94mclose gaps                       \033[0m: (0 if Hybrid assembly, else 1)
        \033[94mmega reads one pass              \033[0m: 0
        \033[94mlimit jump coverage              \033[0m: 300
        \033[94mca parameters                    \033[0m: cgwErrorRate=0.15
        \033[94mjellyfish hash size              \033[0m: based on the estimated size of the genome provided
        \033[94msoap assembly                    \033[0m: 0 (disabled to force CABOG assembly)
        \033[94mflye assembly                    \033[0m: 0 (disabled to force CABOG assembly)
                                    
    \033[92mspades settings
        \033[94mhigh-cov. isolate & multi-cell   \033[0m: --isolate
        \033[94mcoverage cutoff                  \033[0m: auto

    \033[92mflye settings
        \033[94mestimated genome size            \033[0m: est_size
        \033[94mnumber of polishing iterations   \033[0m: 3
        \033[94mcollapse alternative haplotypes  \033[0m: --keep-haplotypes
        
    \033[92mhifasm settings
        \033[0mdefault settings

    \033[96m-----------------------------------------------------------------------\033[0m""")
sleep(250)
log.info("""
    \033[92mracon settings 
        \033[0mdefault settings

    \033[92mpilon settings 
        \033[94mchange file generation           \033[0m: -changes
        \033[94mvcf file generation              \033[0m: -vcf
        \033[94mtracks file generation           \033[0m: -tracks
        \033[94mlargest chunksize limit          \033[0m: 5000000
        \033[94mfix list                         \033[0m: indels, local, snps

    \033[96m-----------------------------------------------------------------------\033[0m""")
sleep(250)
log.info("""
    \033[92mragtag settings 
        \033[0mscaffold
           \033[94mconcatenate unplaced          \033[0m: -C
           \033[94madd suffix to unplaced        \033[0m: -u
        \033[0mcorrect
           \033[94madd suffix to unaltered       \033[0m: -u
        \033[0mpatch
           \033[94madd suffix to unplaced        \033[0m: -u
        
    \033[92mtgs-gapcloser settings
        \033[94mdo not error correct             \033[0m: -ne
            
    \033[92mabyss-sealer settings
        \033[94mpseudoreads length used          \033[0m: 400
        \033[94mbloom filter size                \033[0m: 500M

    \033[96m-----------------------------------------------------------------------\033[0m""")
sleep(250)
log.info("""
    \033[92mfastqc settings
       \033[0mdefault settings
    
    \033[92mnanoplot settings
       \033[94mbivariate plots                   \033[0m: kde, dot
       \033[94mshow logarithmic lengths scaling  \033[0m: --loglength
       \033[94mN50 mark in read length histogram \033[0m: --N50
       \033[94mlog messages to terminal          \033[0m: --verbose

    \033[96m-----------------------------------------------------------------------\033[0m""")
sleep(250)
log.info("""
    \033[92mbusco settings
       \033[94mmode                              \033[0m: genome
       \033[94mforce overwrite                   \033[0m: -f
    
    \033[0m*\033[92mcompleasm settings
       \033[0mdefault settings
    
    \033[92mquast settings
       \033[0mdefault settings

    \033[0m*: Currently disabled since version 3.0.0

\033[91m================================================================================\033[0m
""")
log.info("""
\033[91m.---.\033[92m________\\\033[38;5;208m=\033[96m/\033[92m________\033[91m.---.\033[0m 
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m--------\033[96m/\033[91m=\033[92m\\\033[94m--------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|   \033[94m.\033[96m---------.  \033[94m.\033[96m------.    \033[94m.\033[96m------.    \033[94m.\033[96m-------.\033[0m
\033[91m`---'\033[96m~~~~~~~(\033[38;5;208m===\033[92m)\033[96m~~~~~~~\033[91m`---'  \033[94m/\033[96m|         |\033[94m/\033[96m/        \\ \033[94m/\033[96m/        \\ \033[94m/\033[96m/         \\\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m.--     \033[96m\\\033[91m=\033[92m/    \033[92m,--. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .------\033[94m'\033[96m|   .------\033[94m'\033[96m|   .--\033[94m.   \033[96m. |   .--\033[94m.\033[96m   .\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m|-      \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|  _ \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m----"| \033[96m|   |\033[94m----"| \033[96m|   |\033[94m-'\033[96m|   | |   |\033[94m-'\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[96m`--    \033[92m(\033[91m===\033[96m)   \033[96m`--' \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `----.\033[94m|\033[96m |   |     \033[94m| \033[96m|   +--+   | |   +--+   |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[92m\\\033[38;5;208m=\033[96m/         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|       |\033[94m| \033[96m|   | \033[94m.\033[96m----.|          | |          |\033[0m
 \033[92m|\033[94m|\033[96m|         \033[96m/\033[91m=\033[92m\\         \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  .----'\033[94m| \033[96m|   |\033[94m"\033[96m|    ||   +--+   | |   +------'\033[0m
 \033[92m|\033[94m|\033[96m|        \033[96m(\033[38;5;208m===\033[92m)        \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  |\033[94m---" | \033[96m|   |\033[94m"\033[96m`-.  ||   |\033[94m| \033[96m|   | |   |\033[94m-----"\033[0m
 \033[92m|\033[94m|\033[96m|  \033[96m__     \033[96m\\\033[91m=\033[92m/    \033[96m,__. \033[96m|\033[94m|\033[92m|  \033[94m| \033[96m|  `------.|   `---'  ||   |\033[94m| \033[96m|   | |   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[94m/__\\    \033[92m/\033[38;5;208m=\033[96m\\    \033[94m|__| \033[96m|\033[94m|\033[92m|  \033[94m`.\033[96m|         |`.        \033[96m/\033[94m.\033[96m|   |\033[94m| \033[96m|   |\033[94m.\033[96m|   |\033[0m
 \033[92m|\033[94m|\033[96m| \033[92m|  |   \033[92m(\033[91m===\033[96m)   \033[92m|    \033[96m|\033[94m|\033[92m|    \033[96m`---------'  `------'  `---' \033[94m`\033[96m'---' `---'\033[0m
\033[91m.---.\033[96m_____\033[92m[:\033[94m:\033[96m::::\033[94m:\033[92m]\033[96m_____\033[91m.---.\033[92m    \033[92m╔═══════════════════════════════════════════╗\033[0m
\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|\033[94m------\033[92m|:\033[94m:\033[96m::\033[94m:\033[92m|\033[94m------\033[91m|\033[38;5;208m[\033[93m_\033[38;5;208m]\033[91m|    \033[92m║     \033[94mEnthe\033[96mome Ge\033[97mnome Assemb\033[96mly Pip\033[94meline     \033[92m║\033[0m
\033[91m`---'\033[92m~~~~~~\033[92m|::\033[94m:\033[96m:\033[94m:\033[92m|~~~~~~\033[91m`---'    \033[92m╚═══════════════════════════════════════════╝\033[0m

                    Curated & Maintained by Ian M Bollinger
                         (\033[94mian.bollinger@entheome.org)\033[0m

                              \033[92mdraft_assembly.nf\033[0m
                                version ${version}

\033[91m================================================================================\033[0m
 """)

//
// Processes
//
// Preprocess phase
process preprocess_refseq {
    tag "Preprocess Reference Sequence"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
        
    input:
    val sample_id

    output:
    val true, emit: refseq_preprocess_done
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
        
    echo "DEBUG: sample_id=${sample_id}"
    echo "DEBUG: input_csv=${params.input_csv}"
    echo "DEBUG: output_dir=${params.output_dir}"
    echo "DEBUG: cpu_threads=${params.cpu_threads}"
    echo "DEBUG: ram_gb=${params.ram_gb}"
    
    python3 "${workflow.projectDir}/bin/preprocess_refseq.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process preprocess_illumina {
    tag "Preprocess Illumina Raw Reads"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
        
    input:
    val sample_id
    val refseq_preprocess_done

    output:
    val true, emit: illumina_preprocess_done
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
        
    echo "DEBUG: sample_id=${sample_id}"
    echo "DEBUG: input_csv=${params.input_csv}"
    echo "DEBUG: output_dir=${params.output_dir}"
    echo "DEBUG: cpu_threads=${params.cpu_threads}"
    echo "DEBUG: ram_gb=${params.ram_gb}"
    
    python3 "${workflow.projectDir}/bin/preprocess_illumina.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process preprocess_ont {
    tag "Preprocess ONT Raw Reads"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
    
    input:
    val sample_id
    val illumina_preprocess_done
    
    output:
    val true, emit: ont_preprocess_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    
    echo "DEBUG: sample_id=${sample_id}"
    echo "DEBUG: input_csv=${params.input_csv}"
    echo "DEBUG: output_dir=${params.output_dir}"
    echo "DEBUG: cpu_threads=${params.cpu_threads}"
    echo "DEBUG: ram_gb=${params.ram_gb}"
    
    python3 "${workflow.projectDir}/bin/preprocess_ont.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}"  "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process preprocess_pacbio {
    tag "Preprocess PacBio Raw Reads"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
    
    input:
    val sample_id
    val ont_preprocess_done
        
    output:
    val true, emit: pacbio_preprocess_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    
    echo "DEBUG: sample_id=${sample_id}"
    echo "DEBUG: input_csv=${params.input_csv}"
    echo "DEBUG: output_dir=${params.output_dir}"
    echo "DEBUG: cpu_threads=${params.cpu_threads}"
    echo "DEBUG: ram_gb=${params.ram_gb}"
    
    python3 "${workflow.projectDir}/bin/preprocess_pacbio.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


// Assembly phase
process masurca_assemble {
    tag "MaSuRCA Assembly"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
    
    input:
    val sample_id
    val refseq_preprocess_done
    val illumina_preprocess_done
    val ont_preprocess_done
    val pacbio_preprocess_done

    output:
    val true, emit: masurca_assembly_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    export MPLCONFIGDIR=/tmp/matplotlib-cache
    python3 "${workflow.projectDir}/bin/assemble_masurca.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process spades_assemble {
    tag "SPAdes Assembly"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
    
    input:
    val sample_id
    val masurca_assembly_done

    output:
    val true, emit: spades_assembly_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    export MPLCONFIGDIR=/tmp/matplotlib-cache
    python3 "${workflow.projectDir}/bin/assemble_spades.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process flye_assemble {
    tag "Flye Assembly"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
    
    input:
    val sample_id
    val spades_assembly_done

    output:
    val true, emit: flye_assembly_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    export MPLCONFIGDIR=/tmp/matplotlib-cache
    python3 "${workflow.projectDir}/bin/assemble_flye.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process hifiasm_assemble {
    tag "HiFiasm Assembly"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    maxForks 1
    
    input:
    val sample_id
    val flye_assembly_done

    output:
    val true, emit: hifiasm_assembly_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    export MPLCONFIGDIR=/tmp/matplotlib-cache
    python3 "${workflow.projectDir}/bin/assemble_hifiasm.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process compare_assemblies {
    tag "Compare Assemblies"
    container "${workflow.projectDir}/entheome.sif"
    cpus params.cpu_threads
    memory "${params.ram_gb} GB"
    
    input:
    val sample_id
    val hifiasm_assembly_done

    output:
    val true, emit: compare_assembly_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    python3 "${workflow.projectDir}/bin/compare_assemblies.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process polish_assembly{
    tag "Polish Assembly"
    container "${workflow.projectDir}/entheome.sif"
    maxForks 1
    
    input:
    val sample_id
    val compare_assembly_done
    
    output:
    val true, emit: polish_done
    
    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    python3 "${workflow.projectDir}/bin/polish_assembly.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process curate_assembly{
    tag "Curate Assembly"
    container "${workflow.projectDir}/entheome.sif"
    maxForks 1
    
    input:
    val sample_id
    val polish_done
    
    output:
    val true, emit: curate_done
    
    script:
    """   
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env
    python3 "${workflow.projectDir}/bin/curate_assembly.py" "${sample_id}" "${params.input_csv}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"
    """
}


process final_assembly_qc {
    tag "Final Assembly QC"
    container "${workflow.projectDir}/entheome.sif"
    maxForks 1
    
    input:
    val sample_id
    val curate_done

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh 
    conda activate EGAP_env    
    mkdir -p .temp
    mkdir -p .temp/fontconfig
    export MPLCONFIGDIR=\$PWD/.temp
    export FC_CACHEDIR=\$PWD/.temp/fontconfig

    python3 "${workflow.projectDir}/bin/qc_assessment.py" final "${params.input_csv}" "${sample_id}" "${params.output_dir}" "${params.cpu_threads.toInteger()}" "${params.ram_gb.toInteger()}"    
    """
}


//
// Main workflow
//
workflow {
// Input-Setup
    // Parse the input CSV into sample channels (i.e. one sample_id per row)
    def all_samples_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            [
                row.ONT_SRA             ?: "None",  // [0]
                row.ONT_RAW_DIR         ?: "None",  // [1]
                row.ONT_RAW_READS       ?: "None",  // [2]
                row.ILLUMINA_SRA        ?: "None",  // [3]
                row.ILLUMINA_RAW_DIR    ?: "None",  // [4]
                row.ILLUMINA_RAW_F_READS?: "None",  // [5]
                row.ILLUMINA_RAW_R_READS?: "None",  // [6]
                row.PACBIO_SRA          ?: "None",  // [7]
                row.PACBIO_RAW_DIR      ?: "None",  // [8]
                row.PACBIO_RAW_READS    ?: "None",  // [9]
                row.SAMPLE_ID,                      // [10]
                row.SPECIES_ID,                     // [11]
                row.ORGANISM_KINGDOM,               // [12]
                row.ORGANISM_KARYOTE,               // [13]
                row.BUSCO_1,                        // [14]
                row.BUSCO_2,                        // [15]
                row.EST_SIZE,                       // [16]
                row.REF_SEQ_GCA         ?: "None",  // [17]
                row.REF_SEQ             ?: "None"   // [18]
            ]
        }

// Preprocess
    preprocess_refseq(
        all_samples_ch.map { it[10] } // sample_id
    )

    preprocess_illumina(
        all_samples_ch.map { it[10] }, // sample_id
        preprocess_refseq.out.refseq_preprocess_done
    )

    preprocess_ont(
        all_samples_ch.map { it[10] }, // sample_id
        preprocess_illumina.out.illumina_preprocess_done
    )

    preprocess_pacbio(
        all_samples_ch.map { it[10] }, // sample_id
        preprocess_ont.out.ont_preprocess_done
    )

// Assemble
    // MaSuRCA assemble if Illumina Only or Long Reads (ONT or PacBio) and Illumina are provided
    masurca_assemble(
        all_samples_ch.map { it[10] }, // sample_id
        preprocess_refseq.out.refseq_preprocess_done,
        preprocess_illumina.out.illumina_preprocess_done,
        preprocess_ont.out.ont_preprocess_done,
        preprocess_pacbio.out.pacbio_preprocess_done
    )

    // SPAdes assemble if Illumina Only or Long Reads (ONT or PacBio) and Illumina are provided
    spades_assemble(
        all_samples_ch.map { it[10] }, // sample_id
        masurca_assemble.out.masurca_assembly_done
    )
    
    // Flye assemble if Long Reads (ONT or PacBio) are provided
    flye_assemble(
        all_samples_ch.map { it[10] }, // sample_id
        spades_assemble.out.spades_assembly_done
    )
    
    // HiFiasm assembly if PacBio Reads are provided
    hifiasm_assemble(
        all_samples_ch.map { it[10] }, // sample_id
        flye_assemble.out.flye_assembly_done
    )

// Compare
    compare_assemblies(
        all_samples_ch.map { it[10] }, // sample_id
        hifiasm_assemble.out.hifiasm_assembly_done    
    )

// Polish
    // Polishing waits for Assembly
    polish_assembly(
        all_samples_ch.map { it[10] }, // sample_id
        compare_assemblies.out.compare_assembly_done
    )
    
// Curate
    // Curation waits for Polishish
    curate_assembly(
        all_samples_ch.map { it[10] }, // sample_id
        polish_assembly.out.polish_done
    )

// Assess
    // Assessment waits for Curation
    final_assembly_qc(
        all_samples_ch.map { it[10] }, // sample_id
        curate_assembly.out.curate_done
    )
}
