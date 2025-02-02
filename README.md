# EGAP Pipeline

<div align="center">
  <img src="EGAP_banner.png" alt="EGAP Banner" width="500">
</div>

## Overview

EGAP (Entheome Genome Assembly Pipeline) is a versatile bioinformatics pipeline designed to produce high-quality genome assemblies from **Oxford Nanopore Technologies (ONT)** and/or **Illumina** sequencing data. It can:

1. **Preprocess & QC Reads**  
   - Merge or detect multiple FASTQ files.  
   - Perform read trimming and adapter removal (Trimmomatic, BBDuk).  
   - Deduplicate reads (Clumpify).  
   - Filter and correct ONT reads (Filtlong, Ratatosk).  
   - Generate read metrics (FastQC, NanoPlot, BBMap-based insert-size checks).  

2. **Assembly**  
   - Hybrid or Illumina-only assembly (MaSuRCA), optionally skipping gap closure in MaSuRCA.  
   - If ONT data are available, it also tries Flye or SPAdes, then compares key metrics to pick the best initial assembly.  

3. **Assembly Polishing**
   - Polishes assemblies with Racon (ONT) and Pilon (Illumina).  
   - Optionally removes haplotigs with purge_dups.  

4. **Assembly Curation**
   - Scaffolds and patches using RagTag against a reference genome (if provided).  
   - Performs final gap-closing with either TGS-GapCloser (if ONT) or ABySS-Sealer (if Illumina-only).  

5. **Quality Assessments & Classification**  
   - Runs QUAST for contiguity statistics (N50, L50, GC%, etc.).  
   - Runs Compleasm (BUSCO) on two lineages to measure completeness.  
   - Rates the final assembly as **AMAZING**, **GREAT**, **OK**, or **POOR** based on combined metrics.  
   - (Optional) Integrates coverage calculations against a reference or final assembly size.  

Though optimized for fungal genomes, EGAP can be adapted to many other organisms by switching up lineages, reference sequences, or default thresholds.

**Note on Supported Input Modes:**  
Currently, the pipeline supports only the following combinations:
- **Illumina input only** (SRA, DIR, or RAW FASTQ files)
- **Illumina input + Reference sequence** (either GCA accession or a path to a FASTA file)
- **Illumina input + ONT input** (SRA, DIR, or RAW FASTQ files)
- **Illumina input + ONT input + Reference sequence**

*Future developments include support for:*
- ONT input only
- ONT input + Reference sequence
- PacBio data (SRA, DIR, or RAW FASTQ files) in any combination (solo or with Reference sequence, ONT input, and/or Illumina input)

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Pipeline Flow](#pipeline-flow)
4. [Command-Line Usage](#command-line-usage)
5. [CSV Generation](#csv-generation)
6. [Example Data & Instructions](#example-data--instructions)
7. [Quality Control Output Review](#quality-control-output-review)
8. [Future Improvements](#future-improvements)
9. [References](#references)
10. [Contribution](#contribution)
11. [License](#license)

## Installation

A shell script (e.g., `EGAP_setup.sh`) can install most dependencies (Python 3.8+, Conda, and the main bioinformatics tools):

1. **Run**:

```bash
bash /path/to/EGAP_setup.sh
```

The following tools are installed:
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)
- [BBMap](https://sourceforge.net/projects/bbmap/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [NanoPlot](https://github.com/wdecoster/NanoPlot)
- [Filtlong](https://github.com/rrwick/Filtlong)
- [Ratatosk](https://github.com/DecodeGenetics/Ratatosk)
- [MaSuRCA](https://github.com/alekseyzimin/masurca)
- [Racon](https://github.com/lbcb-sci/racon)
- [Burrows-Wheeler Aligner](https://github.com/lh3/bwa)
- [SamTools](https://github.com/samtools/samtools)
- [BamTools](https://github.com/hartwigmedical/hmftools/tree/master/bam-tools)
- [Pilon](https://github.com/broadinstitute/pilon)
- [purge_dups](https://github.com/dfguan/purge_dups)
- [RagTag](https://github.com/malonge/RagTag)
- [TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser)
- [ABYSS-Sealer](https://github.com/bcgsc/abyss/blob/master/Sealer/sealer.cc)
- [QUAST](https://github.com/ablab/quast)
- [CompleAsm](https://github.com/bioinformatics-centre/compleasm)

Alternatively, you can install the Entheome Ecosystem via Docker. Open a (Linux) terminal in the directory where the `Dockerfile` is located.

```bash
docker build -t entheome_ecosystem .
```

Run the container (adjust the path accordingly):

```bash
docker run -it -v /path/to/data:/path/to/data entheome_ecosystem bash
```

Inside the Docker container, load the pre-generated EGAP environment:

```bash
source /EGAP_env/bin/activate
```

*Future developments include support for:*
- Anaconda installation in a dedicated environment through the Bioconda channel with the following command:

```bash
conda create -y EGAP_env python=3.8 && conda activate EGAP_env && conda install -y -c bioconda egap
```

## Pipeline Flow

<div align="center">
  <img src="EGAP_pipeline.png" alt="EGAP Pipeline" width="500">
</div>

## Command-Line Usage

### Parameters:

- `--input_csv`, `-csv` (str): Path to a CSV containing multiple sample data. (default = None)
- `--ont_sra`, `-osra` (str): Oxford Nanopore Sequence Read Archive (SRA) Accession number. (default = 'SRR27945394')
- `--raw_ont_dir`, `-odir` (str): Path to a directory containing all Raw ONT Reads. (if `-csv` = None; else REQUIRED)
- `--raw_ont_reads`, `-i0` (str): Path to the combined Raw ONT FASTQ reads. (if `-csv` = None; else REQUIRED)
- `--illu_sra`, `-isra` (str): Illumina Sequence Read Archive (SRA) Accession number. (default = 'SRR27945395')
- `--raw_illu_dir`, `-idir` (str): Path to a directory containing all Raw Illumina Reads. (if `-csv` = None; else REQUIRED)
- `--raw_illu_reads_1`, `-i1` (str): Path to the Raw Forward Illumina Reads. (if `-csv` = None; else REQUIRED)
- `--raw_illu_reads_2`, `-i2` (str): Path to the Raw Reverse Illumina Reads. (if `-csv` = None; else REQUIRED)
- `--species_id`, `-ID` (str): Species ID formatted as `<2-letters of Genus>_<full species name>`. (if `-csv` = None; else REQUIRED)
- `--organism_kingdom`, `-Kg` (str): Kingdom the current organism data belongs to. (default: Funga)
- `--organism_karyote`, `-Ka` (str): Karyote type of the organism. (default: Eukaryote)
- `--compleasm_1`, `-c1` (str): Name of the first organism compleasm/BUSCO database to compare to. (default: basidiomycota)
- `--compleasm_2`, `-c2` (str): Name of the second organism compleasm/BUSCO database to compare to. (default: agaricales)
- `--est_size`, `-es` (str): Estimated size of the genome in Mbp (e.g., `60m`). (default: 60m)
- `--ref_seq_gca`, `-rgca` (str): Curated Genome Assembly (GCA) Accession number. (default = None)
- `--ref_seq`, `-rf` (str): Path to the reference genome for assembly. (default: None)
- `--percent_resources`, `-R` (float): Percentage of resources for processing. (default: 1.00)

*Additional PacBio-related parameters (reserved for future support):*
- `--pacbio_sra`, `-psra` (str): PacBio Sequence Read Archive (SRA) Accession number. (default: None)
- `--raw_pacbio_dir`, `-pdir` (str): Path to a directory containing all Raw PacBio Reads. (if `-csv` = None; else REQUIRED)
- `--raw_pacbio_reads`, `-preads` (str): Path to the combined Raw PacBio FASTQ reads. (if `-csv` = None; else REQUIRED)

### Example Command:

```bash
egap --raw_ont_reads /path/to/ont_reads.fq.gz \
     --raw_illu_dir /path/to/illumina_reads/ \
     --species_id AB_speciesname \
     --organism_kingdom Funga \
     --organism_karyote Eukaryote \
     --compleasm_1 basidiomycota \
     --compleasm_2 agaricales \
     --est_size 60m \
     --percent_resources 0.8
```

Or, providing SRA numbers (which will download the files into the current working directory):

```bash
egap --ont_sra SRR######## \
     --illu_sra SRR######## \
     --species_id AB_speciesname \
     --organism_kingdom Funga \
     --organism_karyote Eukaryote \
     --compleasm_1 basidiomycota \
     --compleasm_2 agaricales \
     --est_size 60m \
     --percent_resources 0.8
```

*Note:* Do not use multiple inputs for the same data type (e.g., do NOT use `illu_sra` and `raw_illu_dir` simultaneously).

Alternatively, using a CSV file for multiple samples:

```bash
egap --input_csv /path/to/samples.csv
```

## CSV Generation

To run EGAP with multiple samples, provide a CSV file containing the necessary information for each sample.

### CSV Format

The CSV file should have the following header and columns:

| ONT_SRA       | ONT_RAW_DIR  | ONT_RAW_READS               | ILLUMINA_SRA  | ILLUMINA_RAW_DIR   | ILLUMINA_RAW_F_READS                | ILLUMINA_RAW_R_READS                | SPECIES_ID  | ORGANISM_KINGDOM  | ORGANISM_KARYOTE  | COMPLEASM_1    | COMPLEASM_2  | EST_SIZE  | REF_SEQ_GCA  | REF_SEQ                    |
|---------------|--------------|-----------------------------|---------------|--------------------|-------------------------------------|-------------------------------------|-------------|-------------------|-------------------|----------------|--------------|-----------|--------------|----------------------------|
| None          | None         | /path/to/ONT/sample1.fq.gz  | None          | None               | /path/to/Illumina/sample1_R1.fq.gz  | /path/to/Illumina/sample1_R2.fq.gz  | AB_sample1  | Funga             | Eukaryote         | basidiomycota  | agaricales   | 60m       | None         | None                       |
| None          | None         | None                        | None          | /path/to/Illumina  | None                                | None                                | AB_sample2  | Funga             | Eukaryote         | basidiomycota  | agaricales   | 55m       | None         | /path/to/ref_genome.fasta  |
| SRA2########  | None         | None                        | SRA########   | None               | None                                | None                                | AB_sample3  | Funga             | Eukaryote         | basidiomycota  | agaricales   | 70m       | None         | None                       |

### Column Descriptions

- **ONT_SRA**: Oxford Nanopore Sequence Read Archive (SRA) Accession number. Use `None` if specifying individual files.
- **ONT_RAW_DIR**: Path to the directory containing all Raw ONT Reads. Use `None` if specifying individual files.
- **ONT_RAW_READS**: Path to the combined Raw ONT FASTQ reads (e.g., `/path/to/ONT/sample1.fq.gz`).
- **ILLUMINA_SRA**: Illumina Sequence Read Archive (SRA) Accession number. Use `None` if specifying individual files.
- **ILLUMINA_RAW_DIR**: Path to the directory containing all Raw Illumina Reads. Use `None` if specifying individual files.
- **ILLUMINA_RAW_F_READS**: Path to the Raw Forward Illumina Reads (e.g., `/path/to/Illumina/sample1_R1.fq.gz`).
- **ILLUMINA_RAW_R_READS**: Path to the Raw Reverse Illumina Reads (e.g., `/path/to/Illumina/sample1_R2.fq.gz`).
- **SPECIES_ID**: Species ID formatted as `<2-letters of Genus>_<full species name>` (e.g., `AB_sample1`).
- **ORGANISM_KINGDOM**: Kingdom of the organism (default: `Funga`).
- **ORGANISM_KARYOTE**: Karyote type of the organism (default: `Eukaryote`).
- **COMPLEASM_1**: Name of the first compleasm/BUSCO database (default: basidiomycota).
- **COMPLEASM_2**: Name of the second compleasm/BUSCO database (default: agaricales).
- **EST_SIZE**: Estimated genome size in Mbp (e.g., `60m`).
- **REF_SEQ_GCA**: Curated Genome Assembly (GCA) Accession number (or `None`).
- **REF_SEQ**: Path to the reference genome for assembly (or `None`).

### Example CSV File (`samples.csv`)

```csv
ONT_SRA,ONT_RAW_DIR,ONT_RAW_READS,ILLUMINA_SRA,ILLUMINA_RAW_DIR,ILLUMINA_RAW_F_READS,ILLUMINA_RAW_R_READS,SPECIES_ID,ORGANISM_KINGDOM,ORGANISM_KARYOTE,COMPLEASM_1,COMPLEASM_2,EST_SIZE,REF_SEQ_GCA,REF_SEQ
None,/mnt/d/TESTING_SPACE/Ps_zapotecorum/ONT_MinION/,None,None,/mnt/d/TESTING_SPACE/Ps_zapotecorum/Illumina_PE150/,None,None,Ps_zapotecorum,Funga,Eukaryote,basidiomycota,agaricales,60m,None,None
SRR27945394,None,None,SRR27945395,None,None,None,Ps_caeruleorhiza,Funga,Eukaryote,basidiomycota,agaricales,60m,None,None
None,None,None,SRR5602600,None,None,None,My_speciosa,Flora,Eukaryote,embryophyta,eudicots,200m,GCA_024721245.1,None
```

### Notes

- If you provide a value for `ILLUMINA_RAW_DIR`, set `ILLUMINA_RAW_F_READS` and `ILLUMINA_RAW_R_READS` to `None`. EGAP will automatically detect and process all paired-end reads within that directory. The same applies for `ONT_RAW_DIR`.
- Ensure that all file paths are correct and accessible.
- The CSV file should not contain extra spaces or special characters in the headers.

## Example Data & Instructions

### Create the Folder Structure

First, create the main processing folder with the required sub-folders. Adjust "EGAP_Processing" as needed:

```bash
mkdir -p /path/to/EGAP/EGAP_Processing/ONT /path/to/EGAP/EGAP_Processing/Illumina && cd /path/to/EGAP/EGAP_Processing
```

### Illumina-Only (with Reference Sequence) Assembly Example

*My. speciosa* assembled with reference to its own Reference Sequence.

Download the Reference Sequence into the main processing folder:

```bash
datasets download genome accession GCA_024721245.1 --include genome && unzip ncbi_dataset
```

Download the Illumina data into the Illumina folder (split into multiple files):

```bash
cd Illumina/ && prefetch SRR5602600 && fastq-dump --gzip --split-files SRR5602600 && rm -rf SRR5602600 && cd ..
```

##### Illumina-Only (with Reference Sequence) Assembly Command

Adjust the paths as needed:

```bash
egap --raw_illu_reads_1 /path/to/EGAP/EGAP_Processing/My_speciosa/Illumina/SRR5602600_1.fq.gz \
     --raw_illu_reads_2 /path/to/EGAP/EGAP_Processing/My_speciosa/Illumina/SRR5602600_2.fq.gz \
     --species_id My_speciosa \
     --organism_kingdom Flora \
     --organism_karyote Eukaryote \
     --compleasm_1 embryophyta \
     --compleasm_2 eudicots \
     --est_size 693m \
     --ref_seq /path/to/EGAP/EGAP_Processing/My_speciosa/ncbi_dataset/data/GCA_024721245.1/GCA_024721245.1_ASM2472124v1_genomic.fna
```

### ONT/Illumina Hybrid Assembly Example

*Ps. caeruleorhiza*

Download the ONT data into the ONT folder:

```bash
cd ONT && prefetch SRR13870478 && fastq-dump --gzip SRR27945394 && rm -rf SRR27945394 && cd ..
```

Download the Illumina data into the Illumina folder (split into multiple files):

```bash
cd Illumina && prefetch SRR13870478 && fastq-dump --gzip --split-files SRR27945395 && rm -rf SRR27945395 && cd ..
```

##### ONT/Illumina Hybrid Assembly Command

Adjust the paths as needed:

```bash
egap --raw_ont_reads /path/to/EGAP/EGAP_Processing/ONT/SRR27945394.fastq.gz \
     --raw_illu_reads_1 /path/to/EGAP/EGAP_Processing/Illumina/SRR27945395_1.fastq.gz \
     --raw_illu_reads_2 /path/to/EGAP/EGAP_Processing/Illumina/SRR27945395_2.fastq.gz \
     --species_id Ps_caeruleorhiza \
     --organism_kingdom Funga \
     --organism_karyote Eukaryote \
     --compleasm_1 agaricales \
     --compleasm_2 basidiomycota \
     --est_size 60m
```

## Quality Control Output Review

EGAP generates final assemblies along with:
- **QUAST** metrics (contig count, N50, L50, GC%, coverage)
- **Compleasm (BUSCO)** plots showing Single, Duplicated, Fragmented & Missing scores.
- Final assembly classification: **AMAZING, GREAT, OK,** or **POOR**

### Statistics Thresholds

The current thresholds for each metric classification (subject to change) are:
- **first_compleasm_c** = {"AMAZING": >98.5, "GREAT": >95.0, "OK": >80.0, "POOR": <80.0}
- **second_compleasm_c** = {"AMAZING": >98.5, "GREAT": >95.0, "OK": >80.0, "POOR": <80.0}
- **contigs_thresholds** = {"AMAZING": 100, "GREAT": 1000, "OK": 10000, "POOR": 100000}
- **n50_thresholds** = {"AMAZING": 1000000, "GREAT": 100000, "OK": 1000, "POOR": 100}
- **l50_thresholds** = {"AMAZING": #, "GREAT": #, "OK": #, "POOR": #}

### Compleasm BUSCO Plots

BUSCO outputs are evaluated based on:
- >98.5% Completion (sum of Single and Duplicated genes) for an AMAZING/Great Assembly
- >95.0% Completion for a Good Assembly
- >80% Completion for an OK Assembly
- <80% Completion for a POOR Assembly

Additionally, fewer contigs aligning to BUSCO genes is preferable. Contigs with only duplicated genes are excluded from the plot (noted in the x-axis label).

#### Illumina-Only (with Reference Sequence) Assembly BUSCO Plots

<table align="center">
  <tr>
    <td align="center">
      <img src="My_speciosa_eudicots_odb10_busco.png" alt="My. speciosa eudicots BUSCO plot" width="400">
      <br>
      **My. speciosa eudicots BUSCO**
    </td>
    <td align="center">
      <img src="My_speciosa_embryophyta_odb10_busco.png" alt="My. speciosa embryophyta BUSCO plot" width="400">
      <br>
      **My. speciosa embryophyta BUSCO**
    </td>
  </tr>
</table>

*Flora are known to have large amounts of duplicated genes.*

#### ONT/Illumina Hybrid Assembly BUSCO Plots

<table align="center">
  <tr>
    <td align="center">
      <img src="Ps_caeruleorhiza_agaricales_odb10_busco.png" alt="Ps. caeruleorhiza agaricales BUSCO plot" width="400">
      <br>
      **Ps. caeruleorhiza agaricales BUSCO**
    </td>
    <td align="center">
      <img src="Ps_caeruleorhiza_basidiomycota_odb10_busco.png" alt="Ps. caeruleorhiza basidiomycota BUSCO plot" width="400">
      <br>
      **Ps. caeruleorhiza basidiomycota BUSCO**
    </td>
  </tr>
</table>

*Funga are known to have multinucleate cells and may exhibit higher levels of duplicated genes, though not as extensively as Flora.*

## Future Improvements

- **Automated Quality Assessment Reports**: Generate comprehensive quality reports post-assembly for easier analysis.
- **Improved Data Management**: Automated removal of excess files upon pipeline completion.
- **Enhanced Support for Diverse Genomes**: Optimize pipeline parameters for non-fungal genomes to improve versatility.
- **Improved Error Handling**: Develop robust error detection and user-friendly feedback.
- **Integration with Additional Sequencing Platforms**:  
  - Support for ONT input only and ONT input with Reference sequence.  
  - **Future PacBio Support**: Capability to use PacBio data (SRA, directory, or FASTQ) in any combination (solo or with Reference sequence, ONT input, or Illumina input).

## References

This pipeline was modified from two of the following pipelines:

> Bollinger IM, Singer H, Jacobs J, Tyler M, Scott K, Pauli CS, Miller DR,  
> Barlow C, Rockefeller A, Slot JC, Angel-Mosti V. High-quality draft genomes  
> of ecologically and geographically diverse *Psilocybe* species. *Microbiol Resour Announc* 0:e00250-24; doi: [10.1128/mra.00250-24](https://doi.org/10.1128/mra.00250-24)

> Muñoz-Barrera A, Rubio-Rodríguez LA, Jáspez D, Corrales A, Marcelino-Rodriguez I,  
> Lorenzo-Salazar JM, González-Montelongo R, Flores C. Benchmarking of bioinformatics  
> tools for the hybrid de novo assembly of human whole-genome sequencing data.  
> *bioRxiv* 2024.05.28.595812; doi: [10.1101/2024.05.28.595812](https://doi.org/10.1101/2024.05.28.595812)

The example data are published in:

> Bollinger IM, Singer H, Jacobs J, Tyler M, Scott K, Pauli CS, Miller DR,  
> Barlow C, Rockefeller A, Slot JC, Angel-Mosti V. High-quality draft genomes  
> of ecologically and geographically diverse *Psilocybe* species. *Microbiol Resour Announc* 0:e00250-24; doi: [10.1128/mra.00250-24](https://doi.org/10.1128/mra.00250-24)

> Grassa CJ, Wenger JP, Dabney C, Poplawski SG, Motley ST, Michael TP, Schwartz  
> CJ, Weiblen GD. A complete Cannabis chromosome assembly and adaptive admixture  
> for elevated cannabidiol (CBD) content. *bioRxiv*, 458083; doi: [https://doi.org/10.1101/458083](https://doi.org/10.1101/458083).

## Contribution

If you would like to contribute to the EGAP Pipeline, please submit a pull request or open an issue on GitHub. For major changes, please discuss via an issue first.

## License

This project is licensed under the MIT License.
