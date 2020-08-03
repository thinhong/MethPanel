#!/bin/bash
# Load module R version 3.6.1
module load R/3.6.1
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
methylation_pattern="$BASEDIR/methylation_pattern.r"

outpath1="$OUTPUT/matrix/"
outpath2="$OUTPUT/frequency/"
outpath3="$OUTPUT/pattern/"

mkdir -p "$outpath1" "$outpath2" "$outpath3";
`which R` -f "$methylation_pattern" --args "$INPUT" "$outpath1" "$outpath2" "$outpath3" "$sample" "$AMP" "$FRAME";

echo "1" > "$OUTPUT/${sample}.pattern"

