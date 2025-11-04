## Illumina Mosaic Variant Calling Quick Demo
Here is a quick demo for the Illumina NGS mosaic variant calling using GIAB HG002 chromosome 20 data. The data was sequenced using Illumina NovaSeq system.

```bash
Platform:               Illumina
Sample:     	        GIAB HG002
Input bam coverage:     50x
Control bam coverage:   25x
Reference:              GRCh38_no_alt
Aligner:                BWA-MEM
Region:                 chr20:47000000-47100000
```

**Download data**

```bash
# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/reference_genome/GRCh38_no_alt_chr20.fa
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/reference_genome/GRCh38_no_alt_chr20.fa.fai

# Input and control BAM
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/ilmn/HG002_chr20_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/ilmn/HG002_chr20_demo.bam.bai
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/ilmn/HG003_chr20_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair-mosaic/quick_demo/ilmn/HG003_chr20_demo.bam.bai

INPUT_DIR="${HOME}/ilmn_quick_demo"
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
  --platform ilmn \
  --output_dir ${OUTPUT_DIR} \
  --region chr20:47000000-47100000 \
  --enable_baymgd_tagging \
  --enable_mosaicbase_tagging

### Output
# SNV output: ${INPUT_DIR}/output_paired_mode/snv.vcf.gz
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
  --platform ilmn \
  --output_dir ${OUTPUT_DIR} \
  --region chr20:47000000-47100000 \
  --enable_baymgd_tagging \
  --enable_mosaicbase_tagging
  
### Output
# SNV output: ${INPUT_DIR}/output_single_mode/snv.vcf.gz
```

**Run all commands above:**

```bash
cd ${HOME}
wget "https://raw.githubusercontent.com/HKU-BAL/Clair-Mosaic/main/demo/ilmn_quick_demo.sh"
chmod +x ilmn_quick_demo.sh
./ilmn_quick_demo.sh
```
