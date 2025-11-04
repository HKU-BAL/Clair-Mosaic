<div align="center">
    <img src="images/Clair-Mosaic_icon.png" width="300" alt="Clair-Mosaic">
</div>

# Clair-Mosaic - a deep-learning method for long-read mosaic small variant calling

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Lei Chen, Zhenxian Zheng, Ruibang Luo

Email: lchen@cs.hku.hk, zxzheng@cs.hku.hk, rbluo@cs.hku.hk 

------

## Introduction

Clair-Mosaic is a mosaic variant caller supporting various usage scenarios with both paired samples and single sample as input, and primarily designed for ONT long-read. Clair-Mosaic also applies multiple post-processing strategies, including the usage of germline resources, Bayesian mosaic-germline discriminator, mosaic variant database, and haplotype consistency checking, to distinguish true mosaic variants from germline variants and artifacts.

A preprint describing Clair-Mosaic's designs and results is at [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.10.31.685831v1).

For germline variant calling using **DNA-seq** sample, please try [Clair3](https://github.com/HKU-BAL/Clair3). 

For germline variant calling using **long-read RNA-seq** sample, please try [Clair3-RNA](https://github.com/HKU-BAL/Clair3-RNA).

For somatic variant calling using **tumor-normal pair** samples, please try [ClairS](https://github.com/HKU-BAL/ClairS).

For somatic variant calling using **tumor-only** sample, please try [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO).

------

## Contents
- [Latest Updates](#latest-updates)
- [Quick Demo](#quick-demo)
- [Pre-trained Models](#pre-trained-models)
- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  - [Option 2. Singularity](#option-2-singularity)
  - [Option 3. Docker Dockerfile](#option-3-docker-dockerfile)
- [Usage](#usage)
- [Disclaimer](#disclaimer)

------

## Latest Updates

*v0.0.1 (Sep. 30, 2025)*: Initial release for early access.

------

## Quick Demo

- Oxford Nanopore (ONT) [Q20+](https://nanoporetech.com/q20plus-chemistry) data as input, see [ONT Quick Demo](docs/ont_quick_demo.md).
- PacBio HiFi Revio data as input, see [PacBio HiFi Quick Demo](docs/pacbio_hifi_quick_demo.md).
- Illumina NGS data as input, see [Illumina Quick Demo](docs/illumina_quick_demo.md).

### Quick start

After following [installation](#installation), you can run Clair-Mosaic with one command:

```bash
./run_clair_mosaic --bam_fn input.bam --control_bam_fn control.bam -R ref.fa -o output -t 8 -p ont_r10_dorado_sup_5khz

## Final SNV output VCF file: output/snv.vcf.gz
```

Check [Usage](#usage) for more options.

------

## Pre-trained Models

Clair-Mosaic trained both paired-sample-mode and single-sample-mode models using GIAB samples. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|               Platform               |        Model name         |      Chemistry /Instruments      | Basecaller | Latest update |     Option (`-p/--platform`)      |   Reference   |  Aligner   |
|:------------------------------------:|:-------------------------:|:--------------------------------:|:----------:|:-------------:|:---------------------------------:|:-------------:|:----------:|
|                 ONT                  | r1041_e82_400bps_sup_v420 |          R10.4.1, 5khz           | Dorado SUP | Sep. 30, 2025 |     `ont_r10_dorado_sup_5khz`     | GRCh38_no_alt |  Minimap2  |
|                 ONT                  | r1041_e82_400bps_sup_v410 |          R10.4.1, 4khz           | Dorado SUP | Sep. 30, 2025 |     `ont_r10_dorado_sup_4khz`     | GRCh38_no_alt |  Minimap2  |
|                 ONT                  | r1041_e82_400bps_hac_v410 |          R10.4.1, 4khz           | Dorado HAC | Sep. 30, 2025 |     `ont_r10_dorado_hac_4khz`     | GRCh38_no_alt |  Minimap2  |
|                 ONT                  | r1041_e82_400bps_sup_g615 |          R10.4.1, 4khz           | Guppy6 SUP | Sep. 30, 2025 |     `ont_r10_guppy_sup_4khz`      | GRCh38_no_alt |  Minimap2  |
|               Illumina               |           ilmn            |          NovaSeq/HiseqX          |     -      | Sep. 30, 2025 |              `ilmn`               |    GRCh38     |  BWA-MEM   |
|             PacBio HiFi              |        hifi_revio         | Revio with SMRTbell prep kit 3.0 |     -      | Sep. 30, 2025 |           `hifi_revio`            | GRCh38_no_alt |  Minimap2  |

------

## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available at [DockerHub](https://hub.docker.com/r/hkubal/clair-mosaic). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in docker. 

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair-mosaic:latest \
  /opt/bin/run_clair_mosaic \
    --bam_fn ${INPUT_DIR}/input.bam \            ## use your input bam file name here
    --control_bam_fn ${INPUT_DIR}/control.bam \  ## use your control bam file name here if using paired samples mode; or leave this option empty if using single sample mode
    --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
    --threads ${THREADS} \                       ## maximum threads to be used
    --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_sup_5khz, ont_r10_dorado_sup_4khz, ont_r10_guppy_sup_4khz, ilmn, hifi_revio, etc}
    --output_dir ${OUTPUT_DIR}                   ## output path prefix
```

Check [Usage](#usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in singularity. 

```bash
conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clair-mosaic:latest

# run the sandbox like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clair-mosaic_latest.sif \
  /opt/bin/run_clair_mosaic \
    --bam_fn ${INPUT_DIR}/input.bam \            ## use your input bam file name here
    --control_bam_fn ${INPUT_DIR}/control.bam \  ## use your control bam file name here if using paired samples mode; or leave this option empty if using single sample mode
    --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
    --threads ${THREADS} \                       ## maximum threads to be used
    --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_sup_5khz, ont_r10_dorado_sup_4khz, ont_r10_guppy_sup_4khz, ilmn, hifi_revio, etc}
    --output_dir ${OUTPUT_DIR}                   ## output path prefix
    --conda_prefix /opt/micromamba/envs/clair-mosaic
```

### Option 3. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
git clone https://github.com/HKU-BAL/Clair-Mosaic.git
cd Clair-Mosaic

# build a docker image named hkubal/clair-mosaic:latest
# might require docker authentication to build docker image
docker build -f ./Dockerfile -t hkubal/clair-mosaic:latest .

# run the docker image like option 1
docker run -it hkubal/clair-mosaic:latest /opt/bin/run_clair_mosaic --help
```

------

## Usage

### General Usage

```bash
./run_clair_mosaic \
  --bam_fn ${INPUT_DIR}/input.bam \            ## use your input bam file name here
  --control_bam_fn ${INPUT_DIR}/control.bam \  ## use your control bam file name here if using paired samples mode; or leave this option empty if using single sample mode
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r10_dorado_sup_5khz, ont_r10_dorado_sup_4khz, ont_r10_guppy_sup_4khz, ilmn, hifi_revio, etc}
  --output_dir ${OUTPUT_DIR}                   ## output path prefix
 
## Final SNV output file: ${OUTPUT_DIR}/snv.vcf.gz
```

### Options

**Required parameters:**

```bash
  --bam_fn INPUT_BAM_FN             Input BAM file input. The input file must be samtools indexed.
  -R, --ref_fn FASTA                Reference file input. The input file must be samtools indexed.
  -o, --output_dir OUTPUT_DIR       VCF output directory.
  -t, --threads THREADS             Maximum threads to be used.
  -p, --platform PLATFORM           Select the sequencing platform of the input. Possible options: {ont_r10_dorado_sup_5khz, ont_r10_dorado_sup_4khz, ont_r10_guppy_sup_4khz, ilmn, hifi_revio, etc}.
```

**Miscellaneous parameters:**

```bash
  --control_bam_fn CONTROL_BAM_FN
                        Control BAM file input. The input file must be samtools indexed. Specify the option if using paired samples mode; or leave this option empty if using single sample mode.
  -c CTG_NAME, --ctg_name CTG_NAME
                        The name of the contigs to be processed. Split by ',' for multiple contigs. Default: all contigs will be processed.
  -r REGION, --region REGION
                        A region to be processed. Format: `ctg_name:start-end` (start is 1-based).
  -b BED_FN, --bed_fn BED_FN
                        Path to a BED file. Call variants only in the provided BED regions.
  -G GENOTYPING_MODE_VCF_FN, --genotyping_mode_vcf_fn GENOTYPING_MODE_VCF_FN
                        VCF file input containing candidate sites to be genotyped. Variants will only be called at the sites in the VCF file if provided.
  -H HYBRID_MODE_VCF_FN, --hybrid_mode_vcf_fn HYBRID_MODE_VCF_FN  
                        Enable hybrid calling mode that combines the de novo calling results and genotyping results at the positions in the VCF file given.
  -q QUAL, --qual QUAL  If set, variants with >QUAL will be marked as PASS, or LowQual otherwise.
  --snv_min_af SNV_MIN_AF
                        Minimal SNV AF required for a variant to be called. Decrease SNV_MIN_AF might increase a bit of sensitivity, but in trade of precision, speed and accuracy. Default: 0.05.
  --indel_min_af INDEL_MIN_AF
                        Minimal Indel AF required for a variant to be called. Decrease INDEL_MIN_AF might increase a bit of sensitivity, but in trade of precision, speed and accuracy. Default: 0.05.
  --min_coverage MIN_COVERAGE
                        Minimal coverage required for a variant to be called. Default: 4.
  --chunk_size CHUNK_SIZE
                        The size of each chuck for parallel processing. Default: 5000000.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Define the sample name to be shown in the VCF file. Default: SAMPLE.
  --snv_output_prefix SNV_OUTPUT_PREFIX
                        Prefix for SNV output VCF filename. Default: snv.
  --indel_output_prefix INDEL_OUTPUT_PREFIX
                        Prefix for Indel output VCF filename. Default: indel.
  --remove_intermediate_dir
                        Remove intermediate directory before finishing to save disk space.
  --include_all_ctgs    Call variants on all contigs, otherwise call in chr{1..22} and {1..22}.
  --print_ref_calls     Show reference calls (0/0) in VCF file.
  -d, --dry_run         Print the commands that will be ran.
  --python PYTHON       Absolute path of python, python3 >= 3.9 is required.
  --pypy PYPY           Absolute path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Absolute path of samtools, samtools version >= 1.10 is required.
  --parallel PARALLEL   Absolute path of parallel, parallel >= 20191122 is required.
  --disable_nonmosaic_tagging
                        Disable non-mosaic variants tagging using panel of normals (PoNs). Default: Enabled.
  --enable_baymgd_tagging
                        Enable BayMGD tagging using Bayesian mosaic-germline discriminator. Default: Disabled.
  --enable_mosaicbase_tagging
                        Enable mosaic variant database tagging using MosaicBase. Default: Disabled.
  --enable_post_filtering
                        Enable variants post-processing using haplotype consistency checking. Default: Disabled.
  --enable_indel_calling
                        Enable mosaic Indel calling. Default: Disabled.                   
```

------

## Disclaimer

NOTE: the content of this research code repository (i) is not intended to be a medical device; and (ii) is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.
