#!/bin/bash
mkdir -p fastp_output

for f in ./*.fastq
do
    sample=$(basename "$f" .fastq)

    fastp -i "$f" \
          -o fastp_output/${sample}_trimmed.fastq \
          -q 20 -u 30 \
          --html fastp_output/${sample}_fastp.html \
          --json fastp_output/${sample}_fastp.json \
          --thread 4
done
