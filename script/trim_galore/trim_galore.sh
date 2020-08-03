#!/bin/bash

# Load java version jdk-13.33, python version 3.7.4 and fastqc version 0.11.9
module load java/jdk-13.33
module load python3/3.7.4
module load fastqc/0.11.9
trim_galore="/path/to/trim_galore/0.6.5/trim_galore"

# trimming
if [[ "$isSingleEnd" -eq "1" ]]; then
	# trim single end
	"$trim_galore" -q 30 --illumina --length 20 -o "$O" --clip_R1 1 --three_prime_clip_R1 1 --gzip --fastqc --fastqc_args "--outdir $O" --trim-n "$R1"
	# rename after timmed
	for i in `find "$O" -name "*.fq.gz"`; do mv "$i" "${i/_R1_trimmed.fq.gz/_trimmed_R1.fastq.gz}"; done
else
	# trim paired end
	$trim_galore -q 30 --illumina --length 20 -o "$O" --clip_R1 1 --clip_R2 1 --three_prime_clip_R1 1 --three_prime_clip_R2 1 --gzip --fastqc --fastqc_args "--outdir $O" --trim-n --paired "$R1" "$R2"
	# rename after timmed
	for i in `find "$O" -name "*_R1_trimmed.fq.gz"`; do mv "$i" "${i/_R1_trimmed.fq.gz/_trimmed_R1.fastq.gz}"; done
	for i in `find "$O" -name "*_R2_trimmed.fq.gz"`; do mv "$i" "${i/_R2_trimmed.fq.gz/_trimmed_R2.fastq.gz}"; done
fi

