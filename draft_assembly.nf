#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Check if default input_csv exists, set to null if not. With CLI override capability
params.no_file = params.no_file ?: "$projectDir/assets/NO_FILE" \\ TODO: in processing allow for no_file to occupy input_csv if it is not provided or and _reads paramter not provided or reference_sequence if it is not provided
params.input_csv = params.input_csv ?: "${baseDir}/resources/EGAP_test.csv"
params.input_csv = file(params.input_csv).exists() ? params.input_csv : null

// Input parameters with CLI override capability
params.ont_sra = null
params.raw_ont_dir = null
params.raw_ont_reads = null
params.illu_sra = null
params.raw_illu_dir = null
params.raw_illu_reads_1 = null
params.raw_illu_reads_2 = null
params.pacbio_sra = null
params.raw_pacbio_dir = null
params.raw_pacbio_reads = null
params.species_id = null
params.organism_kingdom = null
params.organism_karyote = null
params.compleasm_1 = null
params.compleasm_2 = null
params.estimated_genome_size = null
params.reference_sequence = null
params.reference_sequence_gca = null
params.percent_resources = 1.0
params.cpu_threads = null
params.ram_gb = null

// Global process resource defaults
process {
    cpus = params.cpu_threads ?: 1
    memory = params.ram_gb ? "${params.ram_gb} GB" : "4 GB"
    container = "${workflow.projectDir}/bin/entheome.sif"
}

log.info("""
.---.________\\=/________.---. 
|[_]|--------/=\\--------|[_]|   .---------.  .------.    .------.    .------.
`---'~~~~~~~(===)~~~~~~~`---'  /|         |//        \\ //        \\ //        \\
 ||| .--     \\=/    ,--. |||  | |  .------'|   .-----/||   .--.   .|   .--.   .
 ||| |-      /=\\    |  _ |||  | |  |     | |   |----' ||   | /|   ||   | /|   |
 ||| `--    (===)   `--' |||  | |  `----.| |   |      ||   +--+   ||   +--+   |
 |||         \\=/         |||  | |       || |   | .----||          ||          |  
 |||         /=\\         |||  | |  .----'| |   |/|    ||   +--+   ||   +------'
 |||        (===)        |||  | |  |     | |   |/`-.  ||   | ||   ||   |-----'
 |||  __     \\=/    ,__. |||  | |  `------.|   `---'  ||   | ||   ||   |
 ||| /__\\    /=\\    |__| |||  \\ |         |`.        /'|   |\\||   ||   |
 ||| |  |   (===)   |    |||    `---------'  `------'  `---'  `---'`---'
.---._____[:::::::]_____.---.  ╔═══════════════════════════════════════════╗
|[_]|------|:::::|------|[_]|  ║     Entheome Genome Assembly Pipeline     ║
`---'~~~~~~|:::::|~~~~~~`---'  ╚═══════════════════════════════════════════╝

              Curated & Maintained by Ian M Bollinger               
                   (ian.bollinger@entheome.org)                     

                          draft_assembly.nf
                                  
FLOWCHART -> 

===============================================================================

    input from   : ${params.base_folder}
    output to    : ${params.output_dir}
    ----------------------------------------------------------------
    run as       : ${workflow.commandLine}
    started at   : ${workflow.start}
    config files : ${workflow.configFiles}
    container    : ${workflow.containerEngine}:${workflow.container}

===============================================================================
""")

// Process to handle input validation and resource calculation
process setup_inputs {
    publishDir "${params.output_dir}/inputs", mode: 'copy'

    input:
    path input_csv from params.input_csv ? Channel.fromPath(params.input_csv, checkIfExists: false) : Channel.fromPath(params.no_file)

    output:
    path "processed_inputs.csv" into input_ch
    val "${file('cpu.txt').text.trim()}" into cpu_ch
    val "${file('ram.txt').text.trim()}" into ram_ch

    script:
    """
    python3 ${baseDir}/process_inputs.py \\
        "${params.input_csv}" \\
        '${params}' \\
        ${params.percent_resources} \\
        "${params.cpu_threads ?: 'None'}" \\
        "${params.ram_gb ?: 'None'}" \\
        ${task.cpus} \\
        ${task.memory.toGiga()}
    """
}

process preprocess_illumina {
    publishDir "${params.output_dir}/preprocess_illumina", mode: 'copy'

    input:
    tuple val(sample), path(csv) from illumina_samples_ch

    output:
    tuple val(sample), path(csv), path("${sample.SPECIES_ID}_*_dedup.fastq.gz") optional true into illumina_preprocessed_ch

    script:
    """
    python3 ${baseDir}/preprocess_illumina.py \\
        "${sample.ILLUMINA_RAW_F_READS}" \\
        "${sample.ILLUMINA_RAW_R_READS}" \\
        "${sample.SPECIES_ID}" \\
        "${sample.RAW_ILLU_DIR}" \\
        "${sample.ILLU_SRA}" \\
        ${task.cpus}
    """
}

process preprocess_ont {
    publishDir "${params.output_dir}/preprocess_ont", mode: 'copy'

    input:
    tuple val(sample), path(csv) from ont_samples_ch
    tuple val(illumina_sample), path(illumina_csv), path(illumina_files) from illumina_preprocessed_ch.ifEmpty([null, null, []])

    output:
    tuple val(sample), path(csv), path(illumina_files), path("${sample.SPECIES_ID}_ont_corrected.fastq.gz") optional true into ont_preprocessed_ch

    script:
    """
    python3 ${baseDir}/preprocess_ont.py \\
        "${sample.ONT_RAW_READS}" \\
        "${sample.EST_SIZE}" \\
        "${sample.SPECIES_ID}" \\
        "${sample.RAW_ONT_DIR}" \\
        "${sample.ONT_SRA}" \\
        "${illumina_files ? illumina_files[0] : 'None'}" \\
        "${illumina_files ? illumina_files[1] : 'None'}" \\
        ${task.cpus}
    """
}

process preprocess_pacbio {
    publishDir "${params.output_dir}/preprocess_pacbio", mode: 'copy'

    input:
    tuple val(sample), path(csv) from pacbio_samples_ch

    output:
    tuple val(sample), path(csv), path([]), path("${sample.SPECIES_ID}_pacbio_filtered.fastq.gz") optional true into pacbio_preprocessed_ch

    script:
    """
    python3 ${baseDir}/preprocess_pacbio.py \\
        "${sample.PACBIO_RAW_READS}" \\
        "${sample.EST_SIZE}" \\
        "${sample.SPECIES_ID}" \\
        "${sample.RAW_PACBIO_DIR}" \\
        "${sample.PACBIO_SRA}" \\
        ${task.cpus}
    """
}

process select_long_reads {
    publishDir "${params.output_dir}/preprocess_long_reads", mode: 'copy'

    input:
    tuple val(sample), path(csv), path(illumina_files), path(ont_file) from ont_preprocessed_ch.ifEmpty([null, null, [], "None"])
    tuple val(sample_pacbio), path(csv_pacbio), path(illumina_files_pacbio), path(pacbio_file) from pacbio_preprocessed_ch.ifEmpty([null, null, [], "None"])
    
    output:
    tuple val(sample), path(csv), path(illumina_files), path("highest_qual_long_reads.fastq.gz") optional true into preprocessed_ch

    when:
    ont_file != "None" || pacbio_file != "None"

    script:
    """
    python3 ${baseDir}/select_long_reads.py \\
        "${ont_file}" \\
        "${pacbio_file}" \\
        "${sample ? sample.SPECIES_ID : sample_pacbio.SPECIES_ID}"
    """
}

process assemble_masurca {
    publishDir "${params.output_dir}/assembly_masurca", mode: 'copy'

    input:
    tuple val(sample), path(csv), path(illumina_files), path(long_reads) from preprocessed_ch
    val ram_gb from setup_inputs.out.ram_ch

    output:
    tuple val(sample), path(csv), path("${sample.SPECIES_ID}_masurca.fasta") optional true into masurca_ch

    when:
    illumina_files.size() > 0  // Requires Illumina reads

    script:
    """
    python3 ${baseDir}/assemble_masurca.py \\
        "${illumina_files[0]}" \\
        "${illumina_files[1]}" \\
        "${long_reads}" \\
        "${sample.SPECIES_ID}" \\
        "${sample.EST_SIZE}" \\
        "${sample.REF_SEQ}" \\
        ${task.cpus} \\
        ${ram_gb}
    """
}

process assemble_spades {
    publishDir "${params.output_dir}/assembly_spades", mode: 'copy'

    input:
    tuple val(sample), path(csv), path(illumina_files), path(long_reads) from preprocessed_ch
    val ram_gb from setup_inputs.out.ram_ch

    output:
    tuple val(sample), path(csv), path("${sample.SPECIES_ID}_spades.fasta") optional true into spades_ch

    when:
    illumina_files.size() > 0  // Requires Illumina reads

    script:
    """
    python3 ${baseDir}/assemble_spades.py \\
        "${illumina_files[0]}" \\
        "${illumina_files[1]}" \\
        "${long_reads}" \\
        "${sample.SPECIES_ID}" \\
        "${sample.REF_SEQ}" \\
        ${task.cpus} \\
        ${ram_gb}
    """
}

process assemble_flye {
    publishDir "${params.output_dir}/assembly_flye", mode: 'copy'

    input:
    tuple val(sample), path(csv), path(illumina_files), path(long_reads) from preprocessed_ch

    output:
    tuple val(sample), path(csv), path("${sample.SPECIES_ID}_flye.fasta") optional true into flye_ch

    when:
    long_reads != "None"  // Requires long reads

    script:
    """
    python3 ${baseDir}/assemble_flye.py \\
        "${long_reads}" \\
        "${sample.SPECIES_ID}" \\
        "${sample.EST_SIZE}" \\
        ${task.cpus}
    """
}

process assemble_hifiasm {
    publishDir "${params.output_dir}/assembly_hifiasm", mode: 'copy'

    input:
    tuple val(sample), path(csv), path(illumina_files), path(long_reads) from preprocessed_ch

    output:
    tuple val(sample), path(csv), path("${sample.SPECIES_ID}_hifiasm.fasta") optional true into hifiasm_ch

    when:
    long_reads != "None" && sample.PACBIO_RAW_READS != "None"  // Requires PacBio reads

    script:
    """
    python3 ${baseDir}/assemble_hifiasm.py \\
        "${long_reads}" \\
        "${sample.SPECIES_ID}" \\
        ${task.cpus}
    """
}

process compare_assemblies {
    publishDir "${params.output_dir}/assembly_comparison", mode: 'copy'

    input:
    tuple val(sample), path(csv), path(illumina_files), path(long_reads) from preprocessed_ch
    tuple val(masurca_sample), path(masurca_csv), path(masurca_assembly) from masurca_ch.ifEmpty([null, null, "None"])
    tuple val(spades_sample), path(spades_csv), path(spades_assembly) from spades_ch.ifEmpty([null, null, "None"])
    tuple val(flye_sample), path(flye_csv), path(flye_assembly) from flye_ch.ifEmpty([null, null, "None"])
    tuple val(hifiasm_sample), path(hifiasm_csv), path(hifiasm_assembly) from hifiasm_ch.ifEmpty([null, null, "None"])

    output:
    tuple val(sample), path(csv), path("${sample.SPECIES_ID}_best_assembly.fasta") into best_assembly_ch

    script:
    """
    python3 ${baseDir}/compare_assemblies.py \\
        "${masurca_assembly}" \\
        "${spades_assembly}" \\
        "${flye_assembly}" \\
        "${hifiasm_assembly}" \\
        "${sample.SPECIES_ID}" \\
        "${sample.COMPLEASM_1}" \\
        "${sample.COMPLEASM_2}" \\
        "${sample.ORGANISM_KARYOTE}" \\
        "${sample.ORGANISM_KINGDOM}" \\
        ${task.cpus}
    """
}

// Updated workflow
workflow {
    setup_inputs()
    setup_inputs.out.input_ch
        .splitCsv(header: true, strip: true)
        .branch {
            illumina: it.ILLUMINA_RAW_F_READS || it.ILLUMINA_RAW_R_READS || it.RAW_ILLU_DIR || it.ILLU_SRA
            ont: it.ONT_RAW_READS || it.RAW_ONT_DIR || it.ONT_SRA
            pacbio: it.PACBIO_RAW_READS || it.RAW_PACBIO_DIR || it.PACBIO_SRA
        }
        .set { branched_samples }

    preprocess_illumina(branched_samples.illumina)
    preprocess_ont(branched_samples.ont, preprocess_illumina.out)
    preprocess_pacbio(branched_samples.pacbio)

    preprocess_ont.out
        .mix(preprocess_pacbio.out)
        .groupTuple(by: 0)
        .map { it -> [it[0], it[1][0], it[2][0], it[3][0] ?: "None", it[3][1] ?: "None"] }
        .set { combined_preprocessed_ch }

    select_long_reads(combined_preprocessed_ch)

    assemble_masurca(select_long_reads.out, setup_inputs.out.ram_ch)
    assemble_spades(select_long_reads.out, setup_inputs.out.ram_ch)
    assemble_flye(select_long_reads.out)
    assemble_hifiasm(select_long_reads.out)

    // Combine assembly outputs for comparison
    select_long_reads.out
        .join(assemble_masurca.out, remainder: true)
        .join(assemble_spades.out, remainder: true)
        .join(assemble_flye.out, remainder: true)
        .join(assemble_hifiasm.out, remainder: true)
        .map { it -> [it[0], it[1], it[2], it[3], it[5] ?: "None", it[7] ?: "None", it[9] ?: "None", it[11] ?: "None"] }
        .set { assembly_comparison_ch }

    compare_assemblies(assembly_comparison_ch)
}