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

# Run Clair-Mosaic in paired-sample mode
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

# Run Clair-Mosaic in single-sample mode
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
