# Concatenate all transcriptome FASTA files into a single genome file
zcat original/REdiscoverTE/rollup_annotation/REdiscoverTE_whole_transcriptome_hg38-20/*.fa.gz > genome.fasta

# Build Salmon index
salmon index -t ./fasta/all_tra.fa --threads 64 -i REdiscoverTE

#!/bin/bash

# Define directories
SAMPLE_DIR="../../data/mRNA_skin/fastq"
OUTPUT_DIR="/home/dell/raid/genome/Ahed/result/gene_salmon/quant_output"  # Output directory

# Create the output directory if it does not exist
mkdir -p $OUTPUT_DIR

# Loop through *_f1.fastq.gz files in the sample directory
for f1 in ${SAMPLE_DIR}/*_f1.fastq.gz; do
  # Get the corresponding *_f2.fastq.gz file
  f2="${f1/_f1.fastq.gz/_f2.fastq.gz}"

  # Ensure that the paired f2 file exists
  if [ -f "$f2" ]; then
    # Extract base filename (without path and extension)
    base=$(basename ${f1%.fastq.gz})  # Remove .fastq.gz suffix
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${base}.salmon.REdiscoverTE"

    # Check if the output directory already exists
    if [ -d "$SAMPLE_OUTPUT_DIR" ]; then
      # Skip processing if the output directory already exists
      echo "Skipping pair: $f1 and $f2, output directory already exists."
    else
      # Process the paired files using Salmon
      echo "Processing pair: $f1 and $f2"
      echo "Base name: $base"

      # Create the sample-specific output directory
      mkdir -p $SAMPLE_OUTPUT_DIR

      # Perform quantification using Salmon for paired-end data
      salmon quant --seqBias --gcBias \
        --index REdiscoverTE \
        --libType A --mates1 ${f1} --mates2 ${f2} \
        --validateMappings \
        -o $SAMPLE_OUTPUT_DIR \
        --threads 64
    fi
  else
    # Warn if no paired f2 file is found
    echo "Warning: No pair found for $f1"
  fi
done

# Generate REdiscoverTE sample metadata file
SAMPLE_DIR=quant_output
echo -e "sample\tquant_sf_path" > ${SAMPLE_DIR}/REdiscoverTE.tsv

ls -1 ${SAMPLE_DIR}/*.REdiscoverTE/quant.sf \
| awk -F/ '{split($(NF-1), a, "_"); print a[1] "\t" $0}'  \
>> ${SAMPLE_DIR}/REdiscoverTE.tsv

# Run REdiscoverTE rollup analysis
REdiscoverTE-master/rollup.R --metadata=./quant_output/REdiscoverTE.tsv \
  --datadir=./data1/ --nozero --threads=64 --assembly=hg38 --outdir=./REdiscoverTE_rollup1/
