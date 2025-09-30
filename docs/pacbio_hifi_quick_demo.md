## PacBio HiFi Mosaic Variant Calling Quick Demo
Here is a quick demo for the PacBio mosaic variant calling using GIAB HG002 chromosome 20 data. The data was sequenced using PacBio latest Revio system.

```bash
Platform:               PacBio HiFi Revio
Sample:     	        GIAB HG002
Input bam coverage:     50x
Control bam coverage:   25x
Reference:              GRCh38_no_alt
Aligner:                pbmm2
Region:                 chr20:47000000-47100000
```

**Download data**

```bash
# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/reference_genome/GRCh38_no_alt_chr20.fa
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/reference_genome/GRCh38_no_alt_chr20.fa.fai

# Input and control BAM
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/pacbio_hifi/HG002_chr20_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/pacbio_hifi/HG002_chr20_demo.bam.bai
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/pacbio_hifi/HG003_chr20_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/pacbio_hifi/HG003_chr20_demo.bam.bai

INPUT_DIR="${HOME}/pacbio_hifi_quick_demo"
mkdir -p ${INPUT_DIR}

REF="GRCh38_no_alt_chr17.fa"
INPUT_BAM="HG002_chr20_demo.bam"
CONTROL_BAM="HG003_chr20_demo.bam"
```

#### Mosaic variant calling using docker pre-built image

##### Run Clair-Mosaic in paired-sample mode

```bash
OUTPUT_DIR="${INPUT_DIR}/output_paired_mode"
mkdir -p ${OUTPUT_DIR}

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair-mosaic:latest \
  /opt/bin/run_clair_mosaic \
  --bam_fn ${INPUT_DIR}/${INPUT_BAM} \
  --control_bam_fn ${INPUT_DIR}/${CONTROL_BAM} \
  --ref_fn ${INPUT_DIR}/${REF} \
  --threads 4 \
  --platform hifi_revio \
  --output_dir ${OUTPUT_DIR} \
  --region chr20:47000000-47100000 \
  --enable_mosaicbase_tagging \
  --enable_post_filtering \
  --enable_indel_calling

### Output
# SNV output: ${INPUT_DIR}/output_paired_mode/snv.vcf.gz
# Indel output: ${INPUT_DIR}/output_paired_mode/indel.vcf.gz
```

##### Run Clair-Mosaic in single-sample mode

```bash
OUTPUT_DIR="${INPUT_DIR}/output_single_mode"
mkdir -p ${OUTPUT_DIR}

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair-mosaic:latest \
  /opt/bin/run_clair_mosaic \
  --bam_fn ${INPUT_DIR}/${INPUT_BAM} \
  --ref_fn ${INPUT_DIR}/${REF} \
  --threads 4 \
  --platform hifi_revio \
  --output_dir ${OUTPUT_DIR} \
  --region chr20:47000000-47100000 \
  --enable_mosaicbase_tagging \
  --enable_post_filtering \
  --enable_indel_calling
  
### Output
# SNV output: ${INPUT_DIR}/output_single_mode/snv.vcf.gz
# Indel output: ${INPUT_DIR}/output_single_mode/indel.vcf.gz
```

**Run all commands above:**

```bash
cd ${HOME}
wget "https://raw.githubusercontent.com/HKU-BAL/Clair-Mosaic/main/demo/pacbio_hifi_quick_demo.sh"
chmod +x pacbio_hifi_quick_demo.sh
./pacbio_hifi_quick_demo.sh
```
