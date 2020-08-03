#!/bin/bash
# Load module R version 3.6.1
module load R/3.6.1
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# project="FieldDefect_14032018_TME"
# path="/g/data3/yo4/Cancer-Epigenetics-Data/Nextseq500/wenjia/$project"
# inpath="$path/polymorphism/"
# outpath="$path/bigTable/"

plot_coverage_methylation="$BASEDIR/plot_coverage_methylation.r"
`which R` -f "$plot_coverage_methylation" --args "$INPATH" "$suffix1" > "$OUTPATH/bigTable${suffix1}.tsv.log"
`which R` -f "$plot_coverage_methylation" --args "$INPATH" "$suffix2" > "$OUTPATH/bigTable${suffix2}.tsv.log"

