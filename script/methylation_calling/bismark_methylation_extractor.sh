#!/bin/bash

# Load modules samtools version 1.9 and bedtools version 2.28.0
module load samtools/1.9
module load bedtools/2.28.0
bismark_methylation_extractor="/path/to/Bismark/0.22.3/bismark_methylation_extractor"

# methylation calling
if [[ "$isSingleEnd" -eq "1" ]]; then
$bismark_methylation_extractor --single-end --output "$OUTPUT" --cytosine_report --comprehensive --merge_non_CpG --genome_folder "$REF" "$BAM";
elif [[ "$isSingleEnd" -eq "0" ]]; then
$bismark_methylation_extractor --paired-end --output "$OUTPUT" --cytosine_report --comprehensive --merge_non_CpG --genome_folder "$REF" "$BAM";	
fi

# make full data frame
infile0="${OUTPUT}/${sample}.bismark.bismark.cov.gz"
infile="${OUTPUT}/${sample}.bismark.cov.gz"
outfile="${OUTPUT}/${sample}.bismark.full.tsv"

mv $infile0 $infile;
echo -e "Amplicon\tCpos\tGpos\tInAmp\tLAmplicon\t${sample}.ratio\t${sample}.C\t${sample}.cov" > "$outfile";
# bedtools intersect -a <(sort -k1,1 -k2,2n -u $frame) -b <(zcat $called| sort -k1,1 -k2,2n -u) -loj| awk '{OFS="\t"}{print $1,$2-200,$3-200,$4,$8,$9,$9+$10}' >> ${called/.bismark.cov.gz/.full.tsv};
for i in `awk '{print $1}' "$FRAME" | sort -u`; do
	if zcat "$infile"| grep -q "$i"; then
		bedtools map -a <(grep "$i" "$FRAME" | sort -k1,1 -k2,2n) -b <(zcat "$infile"| grep "$i"  | sort -k1,1 -k2,2n) -c 4,5,6 -o mean,sum,sum | \
		awk -v add="$add" '{OFS="\t"}{$2=$2-add; $3=$3-add; if($6!="."){$6=100*$7/($7+$8);$8=$7+$8}; print $5,$2,$3,$4,$1,$6,$7,$8}' | \
		sort -k1,1 -k2,2n >> $outfile;
	else
		grep "$i" "$FRAME"| awk '{OFS="\t"}{print $0, ".", ".", "."}' >> "$outfile";
	fi
done