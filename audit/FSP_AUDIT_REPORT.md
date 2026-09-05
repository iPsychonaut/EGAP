# EGAP audit: local vs. GitHub state, and the fspassemblypipeline port proposal

Date: 2026-09-05
Branch: `audit/fsp-port-tests` (local, not pushed; based on `origin/main` @ `3cea58f`)
Reference repo audited: `RBGKew/fspassemblypipeline` @ `4bd5e23` (master, PR #96)

## 1. Local files vs. origin/main

All eleven uploaded root files (`.gitignore`, `build.sh`, `Dockerfile`, `draft_assembly.nf`, `EGAP.py`, `EGAP_setup.sh`, `entheome.sif.def`, `LICENSE`, `meta.yaml`, `nextflow.config`, `README.md`) are content-identical to `origin/main`. `Dockerfile` and `README.md` differ only by CRLF line endings on the local copies. Add a `.gitattributes` with `* text=auto eol=lf` so that stops appearing as a diff.

`bin/` and `resources/` were not uploaded, so their status is unverified. Run `git status --short` locally to close that gap.

Two version inconsistencies exist on main and should be committed as a fix:

`EGAP.py` and `meta.yaml` declare 3.4.2 and the recipe points at `tags/v3.4.2.tar.gz` with a sha256, but the newest remote tag is `v3.4.1`. The bioconda recipe cannot build until `v3.4.2` is tagged on `3cea58f`. The test `test_version_string_matches_recipe` pins EGAP.py and meta.yaml to each other; add the tag check to a release script.

`Dockerfile`, `EGAP_setup.sh`, and `entheome.sif.def` still say 3.4.1. The Dockerfile pulls scripts with `EGAP_BRANCH="v3.4.1"` and an explicit filename list that lacks `preflight_checks.py`, `estimate_runtime.py`, `record_provenance.py`, `monitor_assembly.py`, and `sample_tsv.py`, so a Docker build produces an image whose `EGAP.py` (3.4.2 entry point) imports modules the image does not contain. Fix: `COPY bin/ /EGAP_env/bin/` from the build context instead of a wget list, and derive the version label from `EGAP.py`.

## 2. What fspassemblypipeline actually contains

The email is directionally right about the pipeline's shape and wrong or overstated on four specifics. Each item below was read from source, not from their docs.

### 2.1 Taxonomy-aware BUSCO lineage selection

Location: `subworkflows/local/select_best_assembly_and_qc/main.nf` lines 50 to 81 and 181 to 210 (duplicated verbatim for draft and final runs).

Mechanism: the samplesheet has `taxid,family,order,class,phylum` columns typed by the user. For each sample the code lower-cases each rank name, appends `_odb12`, and returns the first one present in `assets/fungi_busco_lineages.txt` (44 lines, static), else `params.busco_lineage`. The `taxid` column is passed to FCS-GX only; nothing resolves ranks from it. There is no NCBI lookup.

Assessment: worth porting, low cost, but the value is in the whitelist walk, not taxonomy. EGAP's TSV already carries an explicit `BUSCO` pair. The port adds a fallback path for rows where the pair is blank and the rank columns exist. Prototype: `audit/prototypes/busco_lineage.py` (`select_lineage`, `select_lineage_pair`, `whitelist_from_busco_listing`). 10 tests, including both rows from fsp's own reference samplesheet.

Shared library with ANOQI: sensible. The prototype has no imports beyond stdlib and no EGAP dependencies for that reason.

### 2.2 Data-driven k-mer sizing

Location: `modules/local/getseqkitk/main.nf` (awk: `int(len*2/3)`, minus one if even), `modules/local/getkmergeniek/main.nf` plus `bin/extract_best_k.awk` (regex on `<h2>Predicted best k: N</h2>`), `subworkflows/local/genome_assembly/main.nf` lines 43 to 115 (validate odd, 15 <= k <= 127 for SPAdes or 141 for MEGAHIT, not already in defaults, append and sort).

EGAP side: `bin/assemble_spades.py` line 186 hard-codes `["21","33","55","77","99"]` and line 225 hard-codes the matching `K21..K99` cleanup directories. Two copies of one fact.

Assessment: accurate claim. The port is small and the median-read-length heuristic needs no new tool. KmerGenie is an optional extra (bioconda `kmergenie`, R dependency). Prototype: `audit/prototypes/kmer_sizing.py` with a streaming FASTQ median so it does not load the file. 13 tests. One test shows the specific failure the rule prevents: with reads at or under 99 bp the 99-mer must be dropped or SPAdes aborts.

Note: fsp runs three k strategies as three separate assemblies per assembler. For EGAP that triples SPAdes wall time. Recommend one k-list with the predicted k injected, not three branches.

### 2.3 MerquryFK

Location: `select_best_assembly_and_qc/main.nf` runs `FASTK_FASTK` once on cleaned reads, then `MERQURYFK_MERQURYFK` per assembly. Outputs `<prefix>.qv` and `<prefix>.completeness.stats`.

What fsp does not do: use QV for selection. `bin/select_best_assembly.sh` ranks by BUSCO "Complete and single-copy" and breaks ties on QUAST auN. QV is reported to MultiQC only.

EGAP side: `bin/sample_tsv.py` declares `KMER_COMPLETENESS` and `QUAL_VAL` in the stats dict and nothing writes them (confirmed by grep across `bin/`). `compare_assemblies.py` picks by per-metric vote; the tie-break is discussed in section 3.

Assessment: worth doing, and the proposal is stronger than fsp's own use. Tool cost: `fastk` and `merquryfk` from bioconda, plus a FastK k-table per sample (tens of GB for large genomes; delete after use via `file_manager`). Prototype: `audit/prototypes/merqury_qv.py` with header-driven parsers and an explicit lexicographic ranking (`rank_key`) that replaces the vote. 10 tests, including a side-by-side with a faithful copy of EGAP's current vote (`egap_vote`) showing exactly where the two rules disagree.

### 2.4 FCS-GX + Tiara consensus and chimera detection

Location: `subworkflows/local/contamination_detection/main.nf`, `bin/comparison.R` (243 lines), `modules/local/convertfcsrpt/main.nf`.

Mechanism: FCS-GX's `taxonomy.rpt` is cut to a fixed column set, then the R script collapses rows per contig on the id prefix before `~` (FCS-GX itself emits one row per window of a long sequence with a `~start-end` suffix). If windows disagree on species the contig is tagged `possible chimera`; on domain, `possible cross-domain chimera`. Tiara `class_fst_stage` is left-joined afterwards and compared as `match`, `mismatch`, or `chimera`. Contigs under 1 kb are forced to `unknown`. Outputs are a BlobToolKit taxonomy TSV and a synthetic hits TSV.

Nothing is removed. The "consensus" is a labelling step for BlobToolKit. The chimera signal is 100 % FCS-GX; Tiara contributes nothing to it. The R left-join also drops any contig FCS-GX did not report.

Operational cost the email omits: FCS-GX's `gxdb` is roughly 470 GB and NCBI's guidance is to hold it in RAM or on local NVMe. That is a larger footprint than the Kraken2 16 GB standard index EGAP already asks users to provision, and it does not fit the 64 GB "recommended" machine in EGAP's README.

Assessment: the pattern (two classifiers, per-contig agreement table, review file) is worth adopting. The specific second classifier should be a decision, not an assumption. Prototype: `audit/prototypes/contam_consensus.py` reproduces the collapse and compare logic in Python, keeps Tiara-only contigs, and adds the decision layer fsp lacks: remove only on two-classifier agreement, flag on disagreement or chimera, keep otherwise. 30 tests including a 14-row decision matrix. The behavioural change versus EGAP today is that Tiara alone no longer removes a contig; it flags it. Whether that is wanted is a policy call for the maintainers.

Lower-cost alternative for the second classifier: Kraken2 on contigs with the database EGAP already requires. Same consensus code, no new 470 GB dependency.

### 2.5 Selection-script brittleness

`bin/select_best_assembly.sh` derives `reads_type`, `kmer_strategy`, and `assembler` by `cut -d'_'` on filenames, and the QUAST awk block matches on `(R1R2|merged)_`. A sample id containing an underscore, or a future strategy name with two underscores, silently mis-parses. The email's read is correct. EGAP's `compare_assemblies.py` is already a dict of discovered paths; keep it that way.

### 2.6 Checkpoint integrity

EGAP has the same gap the email describes for ANOQI. `bin/decontaminate_assembly.py` line 165 skips Tiara on `os.path.exists(tiara_out)` with no size check, and `write_done_marker` writes free text without verifying that the file it vouches for exists (test `test_done_marker_does_not_verify_outputs`). The same `if os.path.exists(...)` skip pattern appears across the assemblers.

Prototype: `audit/prototypes/checkpoint.py`. `write_done_marker` refuses to write if any declared output is missing or empty, records size (and optional sha256), and writes atomically. `step_can_skip` re-verifies every declared output and rejects legacy free-text markers. 9 tests including the resume-poisoning case (marker written, output truncated afterward). This is the single highest-value item in the email for EGAP because it is a correctness fix, not a feature.

### 2.7 Digest-pinned containers and generated docs

Checked `conf/containers_docker_amd64.config`: four processes (BUSCO, fastp, MerquryFK, MultiQC) use Seqera Wave build-hash tags. Every other module uses nf-core's `quay.io/biocontainers/<tool>:<version>--<build>` tags. Zero `@sha256` digests in the repo. The claim is overstated. The direction (pin, generate docs from one source) is still right and applies to EGAP's own drift in section 1.

### 2.8 Their test discipline

46 nf-test files, 15 for local modules. Spot check: `modules/nf-core/fcsgx/rungx/tests/main.nf.test.snap` asserts `test.taxonomy.rpt:md5,d41d8cd98f00b204e9800998ecf8427e`, the md5 of an empty file. The stub passes by producing nothing. Copy the layout (one test per module, stub mode, CI sharding), not the assertions.

## 3. Defects found in EGAP while writing the baseline tests

`EGAP.find_files` (EGAP.py line 384): `files = []` is rebound by `for root, _, files in os.walk(...)`, and the loop then appends to the list it is iterating. First match hangs forever. No callers anywhere. Delete it. Test runs it in a subprocess with a 5 s timeout and asserts the timeout.

`compare_assemblies` tie-break: `Counter(custom_stats).most_common(1)` where `custom_stats` is built in metric order (BUSCO_1, BUSCO_2, contigs, N50). A 2-2 split therefore goes to the BUSCO_1 winner, not to the first-discovered assembler. Three parametrized cases pin this. The rule is defensible but undocumented; `merqury_qv.rank_key` states it.

`parse_tiara` has no header guard, so `sequence_id -> class_fst_stage` lands in the mapping. Harmless today, wrong by construction.

`get_current_row_data` returns an empty frame for an unknown `SAMPLE_ID`; the failure surfaces later as `IndexError` at `.iloc[0]` inside each stage. Raise at the lookup.

## 4. Recommendation per item

Adopt now, no new tools: checkpoint verification (2.6), explicit ranking rule (2.3 ranking half), k-list derived from read length (2.2), lineage whitelist walk as the fallback when the BUSCO pair is blank (2.1), plus the three small defect fixes in section 3 and the version-drift fix in section 1.

Adopt after a tool-cost decision: MerquryFK QV (adds fastk and merquryfk plus per-sample k-tables). Recommend yes; the QV gap between assemblers is often larger than the BUSCO gap on fungal isolates.

Do not adopt as described: FCS-GX. Keep the consensus code and pick the second classifier separately; Kraken2-on-contigs reuses the database EGAP already needs. Revisit FCS-GX only for users on machines with 512 GB RAM.

Do not adopt: the three-branch k strategy, filename-parsed metadata, and the nf-core template itself. EGAP is a Python orchestrator with a Nextflow wrapper; converting to nf-core is a rewrite, not a port.

## 5. What is in this branch

    .github/workflows/pytest.yml   unit job on py3.8 + py3.11; opt-in conda job gated on RUN_TOOL_TESTS
    pytest.ini, tests/requirements.txt
    tests/conftest.py              sys.path for bin/, mock TSV factory, requires_tools marker
    tests/mockdata.py              seeded generators: FASTA, paired FASTQ, Tiara, FCS-GX rpt,
                                   BUSCO summary, QUAST tsv, MerquryFK qv/stats, KmerGenie html
    tests/test_*.py                35 tests against existing bin/ modules (no bin/ changes)
    tests/prototypes/test_*.py     80 tests against audit/prototypes/
    audit/prototypes/*.py          five candidate modules, stdlib only, not imported by bin/

115 passed on Python 3.8.20 and 3.11.15 in under 8 s. No file under `bin/` or at the repo root was modified.

## 6. Test data

Synthetic data is generated at test time (`tests/mockdata.py`, seeded). For a real end-to-end run without touching SRA, fsp ships two simulated Illumina pairs in `assets/dataset/raw_reads/` (sample_1: 10k pairs at 60x of a small synthetic genome; sample_2: 100k pairs at 40x) and matching pre-built assemblies plus BAMs under `assets/dataset/for_contamination_detection/`. They are small enough to vendor. Alternative public sets: `SRR32496875` (E. coli Illumina, already in `resources/EGAP_test.tsv`) for the shortest real assembly, and BUSCO's `bacteria_odb12` at ~2 MB for a CI-sized lineage.

## 7. Merge gate

Before any prototype moves into `bin/`: the unit job must be green on both Python versions on the PR, the receiving module must gain tests at the same level as the prototype, and one tool-backed run on the fsp `sample_1` data must be attached to the PR with the QC numbers before and after.
