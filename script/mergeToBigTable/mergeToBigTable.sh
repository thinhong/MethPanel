#!/bin/bash

cmd0="paste -d'\t'";
cmd1=$(find "$INPATH" -name "*.bismark.full.tsv"| sort| awk 'NR==1{print "<(cut -f1,2,3,4,5 "$0")"}');
cmd2=$(find "$INPATH" -name "*.bismark.full.tsv"| sort| awk '{print "<(cut -f6,7,8 "$0")"}');
cmd="$cmd0 $cmd1 $cmd2";
eval $cmd > "$OUTPATH/bigTable.tsv"

# filter bigTable
awk '{OFS="\t"}NR==1{print $0}' "$OUTPATH/bigTable.tsv" > "$OUTPATH/bigTable_filtered.tsv"
awk '{OFS="\t"}NR>1{if($4>0) print $0}' "$OUTPATH/bigTable.tsv" >> "$OUTPATH/bigTable_filtered.tsv"

