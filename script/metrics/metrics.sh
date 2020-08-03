#!/bin/bash
# Load module R version 3.6.1
module load R/3.6.1
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

plot_metrics="$BASEDIR/plot_metrics.r"

multiqc_config="$BASEDIR/multiqc.config"
OUTFILE="bigTable/metrics/metrics.tsv"
mkdir -p $(dirname $OUTFILE)

echo -e "Sample\tNumber_Read\tNumber_Mapped_Read\tNumber_Mapped_Converted_Read\tNumber_Mapped_Converted_40mapq_Read\tMethylation_CpG\tMethylation_nonCpG" > $OUTFILE
for j in `find "called" -name "*.bismark_splitting_report.txt"`; do
	i=$(basename "${j/.bismark_splitting_report.txt/}");
	fn1=$(find "aligned/$i" -name "*.bam.1.flagstat"); 
	fn2=$(find "aligned/$i" -name "*.bam.2.flagstat"); 
	fn3=$(find "aligned/$i" -name "*.bam.3.flagstat"); 
	paste <(echo "$i") <(head -n1 "$fn1"| cut -d+ -f1) <(grep " mapped (" "$fn1"| cut -d+ -f1) <(grep " mapped (" "$fn2"| cut -d+ -f1) <(grep " mapped (" "$fn3"| cut -d+ -f1) <(grep "C methylated in CpG context:" "$j"| cut -d: -f2| tr -d '\040\011\012\015') <(grep "C methylated in non-CpG context:" "$j"| cut -d: -f2| tr -d '\040\011\012\015'); 
done| sort -k1,1 >> $OUTFILE

`which R` -f "$plot_metrics" --args "bigTable/metrics/" > "bigTable/metrics/plot_metrics.log"

# trimmed metrics
INPATH="raw_trimmed"
OUTFILE="bigTable/metrics/trimmed_matrix.tsv"
echo -e "Sample\tNo. total reads\tNo. filtered reads (%)\tNo. total basepairs (bp)\t No. Filtered basepairs (%)" > $OUTFILE

for i in $(find $INPATH -name "*_trimming_report.txt"); do 
paste <(echo $i| cut -d"/" -f2) \
<(grep "Total reads processed:" $i| cut -d":" -f2| tr -d -t " "| sed ':a;s/\B[0-9]\{3\}\>/,&/;ta') \
<(grep "Reads written (passing filters):" $i| cut -d"(" -f3| sed "s/%)//g") \
<(grep "Total basepairs processed:" $i| cut -d":" -f2| tr -d -t " "| sed 's/bp//'| sed ':a;s/\B[0-9]\{3\}\>/,&/;ta') \
<(grep "Total written (filtered):" $i| cut -d"(" -f3| sed "s/%)//g"); 
done| sort -k1,1 >> $OUTFILE

# make zip file of all trimmed fastqc for batch downloading  
INPATH="raw_trimmed"
OUTPATH="bigTable/metrics"
mkdir -p "${OUTPATH}/fastqc_html"
find $INPATH -type f -name "*_trimmed_fastqc.html"| xargs -n 1 -I fn cp fn "${OUTPATH}/fastqc_html/"
zip -rm ${OUTPATH}/fastqc_html.zip ${OUTPATH}/fastqc_html

### Accumulating metrics 

# Load python version 3.7.4
module load python3/3.7.4

### metrics for fastq
mkdir -p "bigTable/multiQC/fastq"
cd "bigTable/multiQC/fastq"
# qc report
ls ../../../raw_trimmed/*/*_trimmed_fastqc.zip| xargs -I fn -n 1 unzip fn
# trimming report
ln -s ../../../raw_trimmed/*/*_trimming_report.txt .
# merge qc
multiqc -f -c $multiqc_config .
# remove files
rm -r *_trimmed_fastqc
rm *_trimming_report.txt

### metrics for bam
cd ../../../
mkdir -p "bigTable/multiQC/bam"
cd "bigTable/multiQC/bam"
ln -s ../../../aligned/*/*_SE_report.txt .
ln -s ../../../aligned/*/*.bam.1.flagstat .
ln -s ../../../called/*/*_splitting_report.txt .
ln -s ../../../called/*/*M-bias.txt .
# merge qc
multiqc -f -c $multiqc_config .
# remove files
rm *.txt
rm  *.bam.1.flagstat

