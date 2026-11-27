#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =============================================================================
// EGAP Nextflow runner
// -----------------------------------------------------------------------------
// Runs the same stage scripts, in the same order, as EGAP.py:
//
//   prepare                 EGAP.py --prepare_only: banner, settings (the run's
//                           audit trail, also written to the session log),
//                           pre-flight checks, sample table normalisation
//   preprocess_* ... decontaminate_assembly   one process per bin/<stage>.py
//   qc_assessment, html_reporter, final_compress
//
// tests/test_egap_orchestrator.py fails when this list and EGAP.PROCESSES
// drift apart, so add a stage here whenever one is added to EGAP.py.
//
// Differences from EGAP.py, on purpose:
//   * A failed stage stops the run (Nextflow's default). EGAP.py carries on
//     to QC and the report for that sample. Use -resume after fixing the cause.
//   * Parameters are Nextflow params (nextflow.config or --name value):
//       --input_tsv --output_dir --cpu_threads --ram_gb
//       --dry_run --intermediate_format --busco_only
//       --no_masurca --no_flye --no_spades --no_hifiasm
//
// Usage:
//   nextflow run draft_assembly.nf --input_tsv samples.tsv --output_dir out
// =============================================================================


// Absolute path of a params entry, so every task gets the same stable path.
def absPath(path) {
    return file(path).toAbsolutePath()
}


// Container binds for the sample table's folder and the output folder, so the
// stage scripts can read and write there. Other data drives are bound through
// singularity.runOptions in nextflow.config.
def bindOpts() {
    def tsv_dir = absPath(params.input_tsv).parent
    def out_dir = absPath(params.output_dir)
    return "--bind ${tsv_dir}:${tsv_dir} --bind ${out_dir}:${out_dir}"
}


// Shell preamble shared by every task. Inside the Entheome image EGAP_env has
// to come first on PATH (EGEP_env ships its own python and trimmomatic) and
// CONDA_PREFIX is how preprocess_illumina finds the TruSeq3-PE.fa adapters.
def egapEnv() {
    return """
    set -euo pipefail
    if [ -d /EGAP_env ]; then
        export PATH="/EGAP_env/bin:\$PATH" CONDA_PREFIX=/EGAP_env
    fi
    export MPLCONFIGDIR=/tmp/matplotlib-cache
    export EGAP_INTERMEDIATE_FORMAT="${params.intermediate_format}"
    export EGAP_DRY_RUN="${params.dry_run ? 1 : 0}"
    export EGAP_BUSCO_ONLY="${params.busco_only ? 1 : 0}"
    """.stripIndent()
}


// Command line for one bin/<name>.py stage script (same as EGAP.stage_cmd).
def stage(name, sample_id) {
    def lead = name == 'qc_assessment'
        ? "final \"${absPath(params.input_tsv)}\" \"${sample_id}\""
        : "\"${sample_id}\" \"${absPath(params.input_tsv)}\""
    return egapEnv() + """
    python3 "${workflow.projectDir}/bin/${name}.py" ${lead} \\
        "${absPath(params.output_dir)}" ${params.cpu_threads} ${params.ram_gb}
    """.stripIndent()
}


// True when a sample table cell holds a real value (same as EGAP._cell_has_value).
def hasValue(row, col) {
    def value = row[col]
    return value != null && !(value.toString().trim() in ['', 'None', 'nan'])
}


// QC-only sample: an existing assembly, no reads, no EST_SIZE (same as
// EGAP.sample_is_qc_only). Such a sample only fetches its reference, then
// goes straight to QC and the report.
def isQcOnly(row) {
    def has_ref = ['REF_SEQ', 'REF_SEQ_GCA', 'NT_ASSEMBLY_PATH', 'NT_ASSEMBLY_GCA'].any { col -> hasValue(row, col) }
    def has_reads = row.keySet().any { col ->
        (col.startsWith('ILLUMINA_') || col.startsWith('ONT_') || col.startsWith('PACBIO_')) && hasValue(row, col)
    }
    return has_ref && !has_reads && !hasValue(row, 'EST_SIZE')
}


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
process prepare {
    containerOptions bindOpts()

    output:
    stdout

    script:
    def skips = ['masurca', 'flye', 'spades', 'hifiasm']
        .findAll { name -> params["no_${name}"] }
        .collect { name -> "--no-${name}" }
        .join(' ')
    egapEnv() + """
    python3 "${workflow.projectDir}/EGAP.py" --prepare_only \\
        --input_tsv "${absPath(params.input_tsv)}" --output_dir "${absPath(params.output_dir)}" \\
        --cpu_threads ${params.cpu_threads} --ram_gb ${params.ram_gb} \\
        --intermediate_format ${params.intermediate_format} ${params.dry_run ? '--dry_run' : ''} ${params.busco_only ? '--busco_only' : ''} ${skips}
    """.stripIndent()
}


// -----------------------------------------------------------------------------
// Preprocess
// -----------------------------------------------------------------------------
process preprocess_refseq {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('preprocess_refseq', sample_id)
}

process preprocess_illumina {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('preprocess_illumina', sample_id)
}

process preprocess_ont {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('preprocess_ont', sample_id)
}

process preprocess_pacbio {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('preprocess_pacbio', sample_id)
}

process decontaminate_reads {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('decontaminate_reads', sample_id)
}


// -----------------------------------------------------------------------------
// Assemble (each assembler decides from the sample's read types whether it applies)
// -----------------------------------------------------------------------------
process assemble_masurca {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('assemble_masurca', sample_id)
}

process assemble_flye {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('assemble_flye', sample_id)
}

process assemble_spades {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('assemble_spades', sample_id)
}

process assemble_hifiasm {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('assemble_hifiasm', sample_id)
}


// -----------------------------------------------------------------------------
// Compare, polish, curate, decontaminate
// -----------------------------------------------------------------------------
process compare_assemblies {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('compare_assemblies', sample_id)
}

process polish_assembly {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('polish_assembly', sample_id)
}

process curate_assembly {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('curate_assembly', sample_id)
}

process decontaminate_assembly {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('decontaminate_assembly', sample_id)
}


// -----------------------------------------------------------------------------
// Assess, report, compress
// -----------------------------------------------------------------------------
process qc_assessment {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('qc_assessment', sample_id)
}

process html_reporter {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('html_reporter', sample_id)
}

process final_compress {
    tag "${sample_id}"
    containerOptions bindOpts()

    input:
    val sample_id

    output:
    val sample_id

    script:
    stage('final_compress', sample_id)
}


// -----------------------------------------------------------------------------
// Main workflow
// -----------------------------------------------------------------------------
workflow {
    // One record per sample row. Only SAMPLE_ID travels through the channels;
    // every stage script reads the rest of its row from the table itself.
    def rows = absPath(params.input_tsv).splitCsv(sep: '\t', header: true)
    def sample_ids = rows.collect { row -> row.SAMPLE_ID }
    def qc_only = rows.findAll { row -> isQcOnly(row) }.collect { row -> row.SAMPLE_ID }

    // The output folder has to exist before a container can bind it.
    absPath(params.output_dir).mkdirs()

    // Banner, settings, pre-flight checks and table normalisation run once,
    // before any stage starts.
    def ready = prepare()
    ready.view()

    def refs = preprocess_refseq(ready.flatMap { _settings -> sample_ids }).branch { sample_id ->
        qc: sample_id in qc_only
        build: true
    }

    def ch = refs.build | preprocess_illumina | preprocess_ont | preprocess_pacbio | decontaminate_reads

    // Per-assembler opt-outs, same as EGAP.py --no-<assembler>.
    ch = params.no_masurca ? ch : assemble_masurca(ch)
    ch = params.no_flye    ? ch : assemble_flye(ch)
    ch = params.no_spades  ? ch : assemble_spades(ch)
    ch = params.no_hifiasm ? ch : assemble_hifiasm(ch)

    ch = ch | compare_assemblies | polish_assembly | curate_assembly | decontaminate_assembly

    def reported = ch.mix(refs.qc) | qc_assessment | html_reporter

    // Compress what is left once QC and the report have read it (not on a dry run).
    if( !params.dry_run ) {
        final_compress(reported)
    }
}
