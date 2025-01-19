# EGAP Pipeline

<div align="center">
  <img src="EGAP_banner.png" alt="EGAP Banner" width="500">
</div>

## Overview
EGAP (Entheome Genome Assembly Pipeline) is a versatile bioinformatics pipeline developed for assembling high-quality hybrid genomes using Oxford Nanopore Technologies (ONT) and Illumina sequencing data. It also supports de novo and reference-based assemblies using Illumina data alone. The pipeline encompasses comprehensive steps for read quality control, trimming, genome assembly, polishing, and scaffolding. While optimized for fungal genomes, EGAP can be customized to work with other types of organisms.

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Pipeline Flow](#pipeline-flow)
4. [Command-Line Usage](#command-line-usage)
5. [CSV Generation](#csv-generation)
6. [Example Data & Instructions](#example-data--instructions)
7. [Future Improvements](#future-improvements)
8. [References](#references)

## Installation

The shell script will ensure that Python 3.8 and required libraries are installed. The pipeline has dependencies on a variety of bioinformatics tools, including but not limited to:

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
- [purge_dupes](https://github.com/dfguan/purge_dups)
- [RagTag](https://github.com/malonge/RagTag)
- [TGS-GapCloser](https://github.com/BGI-Qingdao/TGS-GapCloser)
- [ABYSS-Sealer](https://github.com/bcgsc/abyss/blob/master/Sealer/sealer.cc)
- [QUAST](https://github.com/ablab/quast)
- [CompleAsm](https://github.com/bioinformatics-centre/compleasm)
- [Merqury](https://github.com/marbl/merqury)

You can install pre-requisites using the shell script:

```bash
bash /path/to/EGAP_setup.sh
```
This creates a conda environment named "EGAP_env".

Alternatively you can install the Entheome Ecosystem via Docker. Open a (Linux) terminal in the directory where the "Dockerfile" is located.
```bash
docker build -t entheome_ecosystem .
```

Run the Container, change "path/to/data" to match appropriate pathing.
```bash
docker run -it -v /path/to/data:/path/to/data entheome_ecosystem bash
```

In the Docker Image, load the pre-generated EGAP_env.
```bash
source /EGAP_env/bin/activate
```


## Pipeline Flow

<div align="center">
  <img src="EGAP_pipeline.png" alt="EGAP Pipline" width="500">
</div>

## Command-Line Usage

### Parameters:

- \`--input_csv\`, \`-csv\` (str): Path to a CSV containing multiple sample data. (default = None)
- \`--raw_ont_dir\`, \`-odir\` (str): PPath to a directory containing all Raw ONT Reads. (if \`-csv\` = None; else REQUIRED)
- \`--raw_ont_reads\`, \`-i0\` (str): Path to the combined Raw ONT FASTQ reads. (if \`-csv\` = None; else REQUIRED)
- \`--raw_illu_dir\`, \`-idir\` (str): Path to a directory containing all Raw Illumina Reads. (if \`-csv\` = None; else REQUIRED)
- \`--raw_illu_reads_1\`, \`-i1\` (str): Path to the Raw Forward Illumina Reads. (if \`-csv\` = None; else REQUIRED)
- \`--raw_illu_reads_2\`, \`-i2\` (str): Path to the Raw Reverse Illumina Reads. (if \`-csv\` = None; else REQUIRED)
- \`--species_id\`, \`-ID\` (str): Species ID formatted as \`<2-letters of Genus>_<full species name>\`. (if \`-csv\` = None; else REQUIRED)
- \`--organism_kingdom\`, \`-Kg\` (str): Kingdom the current organism data belongs to. (default: Funga)
- \`--organism_karyote\`, \`-Ka\` (str): Karyote type of the organism. (default: Eukaryote)
- \`--compleasm_1\`, \`-c1` (str): Name of the first organism compleasm/BUSCO database to compare to. (default: basidiomycota)
- \`--compleasm_2\`, \`-c2` (str): Name of the second organism compleasm/BUSCO database to compare to. (default: agaricales)
- \`--est_size\`, \`-es\` (str): Estimated size of the genome in Mbp (million base pairs). (default: 60m)
- \`--ref_seq\`, \`-rf\` (str): Path to the reference genome for assembly. (default: None)
- \`--percent_resources\`, \`-R\` (float): Percentage of resources for processing. (default: 1.00)

### Example Command:

```bash
python /path/to/EGAP.py --raw_ont_reads /path/to/ont_reads.fq.gz \
                        --raw_illu_dir /path/to/illumina_reads/ \
                        --species_id AB_speciesname \
                        --organism_kingdom Funga \
                        --organism_karyote Eukaryote \
                        --compleasm_1 basidiomycota \
                        --compleasm_2 agaricales \
                        --est_size 60m \
                        --percent_resources 0.8
```

Alternatively, using a CSV file for multiple samples:

```bash
python /path/to/EGAP.py --input_csv /path/to/samples.csv
```

## CSV Generation

To run EGAP with multiple samples, you can provide a CSV file containing the necessary information for each sample. Below is the correct format for the CSV file:

### CSV Format

The CSV file should have the following header and columns:

| ONT_RAW_DIR   | ONT_RAW_READS                  | ILLUMINA_RAW_DIR   | ILLUMINA_RAW_F_READS                | ILLUMINA_RAW_R_READS                | SPECIES_ID     | ORGANISM_KINGDOM | ORGANISM_KARYOTE | COMPLEASM_1   | COMPLEASM_2 | EST_SIZE | REF_SEQ                    |
|---------------|--------------------------------|--------------------|-------------------------------------|-------------------------------------|----------------|------------------|------------------|---------------|-------------|----------|----------------------------|
| None          | /path/to/ONT/sample1.fq.gz     | None               | /path/to/Illumina/sample1_R1.fq.gz  | /path/to/Illumina/sample1_R2.fq.gz  | AB_sample1     | Funga            | Eukaryote        | basidiomycota | agaricales  | 60m      | /path/to/ref_genome1.fasta |
| /path/to/ONT  | None                           | /path/to/Illumina  | None                                | None                                | AB_sample2     | Funga            | Eukaryote        | basidiomycota | agaricales  | 55m      | /path/to/ref_genome2.fasta |

### Column Descriptions

- **ONT_RAW_DIR**: Path to the directory containing all Raw ONT Reads. Use `None` if specifying individual read files.
- **ONT_RAW_READS**: Path to the combined Raw ONT FASTQ reads (e.g., `/path/to/ONT/sample1.fq.gz`).
- **ILLUMINA_RAW_DIR**: Path to the directory containing all Raw Illumina Reads. Use `None` if specifying individual read files.
- **ILLUMINA_RAW_F_READS**: Path to the Raw Forward Illumina Reads (e.g., `/path/to/Illumina/sample1_R1.fq.gz`).
- **ILLUMINA_RAW_R_READS**: Path to the Raw Reverse Illumina Reads (e.g., `/path/to/Illumina/sample1_R2.fq.gz`).
- **SPECIES_ID**: Species ID formatted as `<2-letters of Genus>_<full species name>` (e.g., `AB_sample1`).
- **ORGANISM_KINGDOM**: Kingdom the current organism data belongs to (default: `Funga`).
- **ORGANISM_KARYOTE**: Karyote type of the organism. (default: Eukaryote).
- **COMPLEASM_1**: Name of the first organism compleasm/BUSCO database to compare to. (default: basidiomycota).
- **COMPLEASM_2**: Name of the second organism compleasm/BUSCO database to compare to. (default: agaricales).
- **EST_SIZE**: Estimated size of the genome in Mbp (million base pairs) (e.g., `60m`).
- **REF_SEQ**: Path to the reference genome for assembly. Use `None` if not applicable.

### Example CSV File (`samples.csv`)

```csv
ONT_RAW_DIR,ONT_RAW_READS,ILLUMINA_RAW_DIR,ILLUMINA_RAW_F_READS,ILLUMINA_RAW_R_READS,SPECIES_ID,ORGANISM_KINGDOM,ORGANISM_KARYOTE,COMPLEASM_1,COMPLEASM_2,EST_SIZE,REF_SEQ
None,/mnt/d/EGAP/EGPA_Processing/Ps_zapotecorum/ONT/SRR########.fastq.gz,None,/mnt/d/EGAP/EGAP_Processing/Ps_zapotecorum/Illumina/SRR########_1.fq.gz,/mnt/d/EGAP/EGAP_Processing/Ps_zapotecorum/Illumina/SRR########_2.fq.gz,Ps_zapotecorum,Funga,Eukaryote,basidiomycota,agaricales,60m,None
None,/mnt/d/EGAP/EGAP_Processing/Ps_gandalfiana/ONT/SRR########.fastq.gz,/mnt/d/EGAP/EGAP_Processing/Ps_gandalfiana/Illumina/B1_3,None,None,Ps_gandalfiana,Funga,Eukaryote,basidiomycota,agaricales,60m,/mnt/d/EGAP/EGAP_Processing/Ps_gandalfiana/GCF_#########_#.fna
```

### Notes

- If you provide a value for `ILLUMINA_RAW_DIR`, set `ILLUMINA_RAW_F_READS` and `ILLUMINA_RAW_R_READS` to `None`. EGAP will automatically detect and process all paired-end reads within the specified directory. This is also True if for if you provide `ONT_RAW_DIR`.
- Ensure that all file paths are correct and accessible.
- The CSV file should not contain any extra spaces or special characters in the headers.

## Example Data & Instructions

### Create the Folder Structure

First, create the main processing folder with the required sub-folders; change "EGAP_Processing" or add your own organism specific folder as needed:

```bash
mkdir -p /path/to/EGAP/EGAP_Processing/ONT /path/to/EGAP/EGAP_Processing/Illumina && \
cd /path/to/EGAP/EGAP_Processing
```

### Illumina-Only (with Reference Sequence) Assembly Example:

My. speciosa assembled with reference to the same organism's Reference Sequence.

Download the Reference Sequence into main processing folder:

```bash
datasets download genome accession GCA_024721245.1 --include genome && \
unzip ncbi_dataset
```

Download the Illumina data into the Illumina folder (split into multiple files):

```bash
cd Illumina/ && \
prefetch SRR5602600 && \
fastq-dump --gzip --split-files SRR5602600 && \
rm -rf SRR5602600 && \
cd ..
```

##### Illumina-Only (with Reference Sequence) Assembly Command:

Adjust the paths to correctly match the downloaded files.

```bash
python /mnt/d/EGAP/EGAP.py --raw_illu_reads_1 /path/to/EGAP/EGAP_Processing/My_speciosa/Illumina/SRR5602600_1.fq.gz \
                           --raw_illu_reads_2 /path/to/EGAP/EGAP_Processing/My_speciosa/Illumina/SRR5602600_2.fq.gz \
                           --species_id My_speciosa \
                           --organism_kingdom Flora \
                           --organism_karyote Eukaryote \
                           --compleasm_1 embryophyta \
                           --compleasm_2 eudicots \
                           --est_size 693m \
                           --ref_seq /path/to/EGAP/EGAP_Processing/My_speciosa/ncbi_dataset/data/GCA_024721245.1/GCA_024721245.1_ASM2472124v1_genomic.fna
```


### ONT/Illumina Hybrid Assembly Example: 

Ps. caeruleorhiza

Download the ONT data into the ONT folder:

```bash
cd ONT && \
prefetch SRR13870478 && \
fastq-dump --gzip SRR27945394 && \
rm -rf SRR27945394 && \
cd ..
```

Download the Illumina data into the Illumina folder (split into multiple files):

```bash
cd Illumina && \
prefetch SRR13870478 && \
fastq-dump --gzip --split-files SRR27945395 && \
rm -rf SRR27945395 && \
cd ..
```

#####  ONT/Illumina Hybrid Assembly Command:

Adjust the paths to correctly match the downloaded files.

```bash
python /mnt/d/EGAP/EGAP.py --raw_ont_reads /path/to/EGAP/EGAP_Processing/ONT/SRR27945394.fastq.gz \
                           --raw_illu_reads_1 /path/to/EGAP/EGAP_Processing/Illumina/SRR27945395_1.fastq.gz \
                           --raw_illu_reads_2 /path/to/EGAP/EGAP_Processing/Illumina/SRR27945395_2.fastq.gz \
                           --species_id Ps_caeruleorhiza \
                           --organism_kingdom Funga \
                           --organism_karyote Eukaryote \
                           --compleasm_2 basidiomycota \
                           --compleasm_1 agaricales \
                           --est_size 60m
```

## Quality Control Output Review
EGAP Will rate the final assembly based on QUAST & BUSCO statistics. The final assessment will be one of the following:
- **AMAZING**
- **GREAT**
- **OK**
- **POOR**

### QUAST Statistics
(TBD)
- **contigs_thresholds** = {"AMAZING": 100, "GREAT": 1000, "OK": 10000, "POOR": 100000}
- **n50_thresholds** = {"AMAZING": 1000000, "GREAT": 100000, "OK": 1000, "POOR": 100}
- **l50_thresholds** = {"AMAZING": #, "GREAT": #, "OK": #, "POOR": #}

### Compleasm BUSCO Plots 
When assessing BUSCO outputs we desire to see >98.5% Completion (Sum of Single and Duplicated genes) in an AMAZING GREAT Assembly; >95% Completion in a Good Assembly; >90% Completion in an OK Assembly; and <90% Completion in a POOR Assembly. Along with High Completion, it is desired to see very few contigs that the genes align to; these sequences can be indiciative of chromsome-candidates. Any sequences that only matched Duplicated genes were excluded from the plot (but noted in the x-axis label).

#### Illumina-Only (with Reference Sequence) Assembly BUSCO Plots:
<table align="center">
  <tr>
    <td align="center">
      <img src="My_speciosa_eudicots_odb10_busco.png" alt="My. speciosa eudicots BUSCO plot" width="400">
      <br>
      <b>My. speciosa eudicots BUSCO</b>
    </td>
    <td align="center">
      <img src="My_speciosa_embryophyta_odb10_busco.png" alt="My. speciosa embryophyta BUSCO plot" width="400">
      <br>
      <b>My. speciosa embryophyta BUSCO</b>
    </td>
  </tr>
</table>
Flora are known to have large amounts of duplicated genes.

#### ONT/Illumina Hybrid Assembly BUSCO Plots:
<table align="center">
  <tr>
    <td align="center">
      <img src="Ps_caeruleorhiza_agaricales_odb10_busco.png" alt="Ps. caeruleorhiza agaricales BUSCO plot" width="400">
      <br>
      <b>Ps. caeruleorhiza agaricales BUSCO</b>
    </td>
    <td align="center">
      <img src="Ps_caeruleorhiza_basidiomycota_odb10_busco.png" alt="Ps. caeruleorhiza basidiomycota BUSCO plot" width="400">
      <br>
      <b>Ps. caeruleorhiza basidiomycota BUSCO</b>
    </td>
  </tr>
</table>
Funga are known to have multi-nucleate cells, and thus have a higher possibility of duplicated genes, but not to the same degree as Flora. 


## Future Improvements
- **Automated Quality Assessment Reports**: Generate comprehensive quality reports post-assembly for easier analysis.
- **Improved Data Management**: Removal of excess files once pipeline complete.
- **Enhanced Support for Diverse Genomes**: Optimize pipeline parameters for non-fungal genomes to improve versatility.
- **Improved Error Handling**: Develop more robust error detection and user-friendly feedback mechanisms.
- **Integration with Additional Sequencing Platforms**: Expand support beyond ONT and Illumina to include platforms like PacBio.

## References
This pipeline was modified From two of the following pipelines:

    Bollinger IM, Singer H, Jacobs J, Tyler M, Scott K, Pauli CS, Miller DR,
    Barlow C, Rockefeller A, Slot JC, Angel-Mosti V. High-quality draft genomes
    of ecologically and geographically diverse Psilocybe species. Microbiol Resour
    Announc 0:e00250-24; doi: https://doi.org/10.1128/mra.00250-24
    
    Muñoz-Barrera A, Rubio-Rodríguez LA, Jáspez D, Corrales A , Marcelino-Rodriguez I,
    Lorenzo-Salazar JM, González-Montelongo R, Flores C. Benchmarking of bioinformatics
    tools for the hybrid de novo assembly of human whole-genome sequencing data.
    bioRxiv 2024.05.28.595812; doi: https://doi.org/10.1101/2024.05.28.595812 

The example data are published in:

    Bollinger IM, Singer H, Jacobs J, Tyler M, Scott K, Pauli CS, Miller DR,
    Barlow C, Rockefeller A, Slot JC, Angel-Mosti V. High-quality draft genomes
    of ecologically and geographically diverse Psilocybe species. Microbiol Resour
    Announc 0:e00250-24; doi: https://doi.org/10.1128/mra.00250-24

    McKernan K, Kane L, Helbert Y, Zhang L, Houde N, McLaughlin S. A whole genome
    atlas of 81 Psilocybe genomes as a resource for psilocybin production. F1000Research
    2021, 10:961; doi: https://doi.org/10.12688/f1000research.55301.2
    
    Grassa CJ, Wenger JP, Dabney C, Poplawski SG, Motley ST, Michael TP, Schwartz
    CJ, Weiblen GD. A complete Cannabis chromosome assembly and adaptive admixture
    for elevated cannabidiol (CBD) content. BioRxiv, 458083; doi:https://doi.org/10.1101/458083.

## Contribution
If you would like to contribute to the EGAP Pipeline, please submit a pull request or open an issue on GitHub. For major changes, please discuss them with us first via an issue.

## License
This project is licensed under the MIT License.
