#!/bin/bash
# Load R 3.6.1
module load R/3.6.1
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# project="FieldDefect_14032018_TME"
# inpath="/g/data3/yo4/Cancer-Epigenetics-Data/Nextseq500/wenjia/$project/pattern/pattern"
# AMP="/g/data3/yo4/Cancer-Epigenetics-Data/Nextseq500/wenjia/$project/config/TME_pan12_name_sorted_mapped.bed"
# outpath="/g/data3/yo4/Cancer-Epigenetics-Data/Nextseq500/wenjia/$project/polymorphism"
# mkdir -p "$outpath"

mergeToAmplicon="$BASEDIR/mergeToAmplicon.r"

# i="PC51_FERD3L::chr7:19184732-19185193"
for i in `awk '{print $8"::"$5":"$6"-"$7}' "$AMP"`; do
	echo "$i"
	`which R` -f "$mergeToAmplicon" --args "$INPATH" "$i" "$OUTPATH";
done &> "$OUTPATH/Rout.txt"

