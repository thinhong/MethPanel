#!/bin/bash

# Load bowtie2 version 2.3.5.1 and samtools version 1.9
module load bowtie2/2.3.5.1
module load samtools/1.9

bismark="/path/to/Bismark/0.22.3/bismark"
filter="/path/to/Bismark/0.22.3/filter_non_conversion"
R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
mkdir -p "$OUTPUT/temp" 

# -B Write all output to files starting with this base file name
# $bismark/bismark --bowtie2 --non_directional -o $3 -N 1 $2 -1 ${1}_R1.fastq.gz -2 ${1}_R2.fastq.gz

if [[ "$isSingleEnd" -eq "1" && "$REMOVE_UNCONVERTED" -eq "1" ]] ; 
	then
		$bismark --bowtie2 --non_directional --temp_dir "$OUTPUT/temp" -o "$OUTPUT" -N 1 -l 28 -B "${sample}" "$REF" "$R1";
		samtools flagstat "$OUTPUT/${sample}.bam" > "$OUTPUT/${sample}.bam.1.flagstat";
		$filter -s --threshold "$NUMBER_C" "$OUTPUT/${sample}.bam";
		sleep 5;
		samtools view -q 40 -h "$OUTPUT/${sample}.nonCG_filtered.bam" | samtools sort -T "$OUTPUT/temp" -@8 - > "$OUTPUT/${sample}.bismark.bam";
		samtools flagstat "$OUTPUT/${sample}.nonCG_filtered.bam" > "$OUTPUT/${sample}.nonCG_filtered.bam.2.flagstat";
elif [[ "$isSingleEnd" -eq "1" && "$REMOVE_UNCONVERTED" -eq "0" ]] ; 
	then
		$bismark --bowtie2 --non_directional --temp_dir "$OUTPUT/temp" -o "$OUTPUT" -N 1 -l 28 -B "${sample}" "$REF" "$R1";
		sleep 5;
		samtools flagstat "$OUTPUT/${sample}.bam" > "$OUTPUT/${sample}.bam.1.flagstat";
		samtools view -q 40 -h "$OUTPUT/${sample}.bam" | samtools sort -T "$OUTPUT/temp" -@8 - > "$OUTPUT/${sample}.bismark.bam";
elif [[ "$isSingleEnd" -eq "0" && "$REMOVE_UNCONVERTED" -eq "1" ]] ; 
	then
		$bismark --bowtie2 --non_directional --temp_dir "$OUTPUT/temp" -o "$OUTPUT" -N 1 -l 28 -B "${sample}" "$REF" -1 "$R1" -2 "$R2";
		samtools flagstat "$OUTPUT/${sample}.bam" > "$OUTPUT/${sample}.bam.1.flagstat";
        $filter -p --threshold "$NUMBER_C" "$OUTPUT/${sample}.bam";
        sleep 5;
		samtools view -q 40 -h "$OUTPUT/${sample}.nonCG_filtered.bam" | samtools sort -T "$OUTPUT/temp" -@8 - > "$OUTPUT/${sample}.bismark.bam";
		samtools flagstat "$OUTPUT/${sample}.nonCG_filtered.bam" > "$OUTPUT/${sample}.nonCG_filtered.bam.2.flagstat";
elif [[ "$isSingleEnd" -eq "0" && "$REMOVE_UNCONVERTED" -eq "0" ]] ; 
	then
		$bismark --bowtie2 --non_directional --temp_dir "$OUTPUT/temp" -o "$OUTPUT" -N 1 -l 28 -B "${sample}" "$REF" -1 "$R1" -2 "$R2";
		sleep 5;
		samtools flagstat "$OUTPUT/${sample}.bam" > "$OUTPUT/${sample}.bam.1.flagstat";
		samtools view -q 40 -h "$OUTPUT/${sample}.bam" | samtools sort -T "$OUTPUT/temp" -@8 - > "$OUTPUT/${sample}.bismark.bam";
fi

# index
samtools index "$OUTPUT/${sample}.bismark.bam";
# flagstat
samtools flagstat "$OUTPUT/${sample}.bismark.bam" > "$OUTPUT/${sample}.bismark.bam.3.flagstat";
# clean the temp
rm -r "$OUTPUT/temp"

