#!/bin/bash
# Load module R version 3.6.1
module load R/3.6.1
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

mergeToBigTableFullRead="$BASEDIR/mergeToBigTableFullRead.r"
`which R` -f "$mergeToBigTableFullRead" --args "$INPATH" "$OUTPATH" > "$OUTPATH/bigTable_fullRead.tsv.mergeToBigTableFullRead.log"

