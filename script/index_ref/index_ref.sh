####### Amplicon regions

# Load bowtie2 version 2.3.5.1, bedtools version 2.28.0 and R version 3.6.1
module load bowtie2/2.3.5.1
module load bedtools/2.28.0
module load R/3.6.1
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bismark="/path/to/Bismark/0.22.3"
create_CpG_ref="$BASEDIR/create_CpG_ref.r"
mkdir -p $REF

# name the fasta
build_name_panel="${build_name}_${panel_name}"
refs="$REF/${build_name_panel}.fa"
# name the amplicon coordinate with .extend.merged
inAmp_extend_merged_bed="$REF/${panel_name}.extend.merged.bed"
# name the amplicon coordinate with _sorted_mapped
inAmp_sorted_mapped="$REF/${panel_name}_sorted_mapped.bed"

# project="FieldDefect_14032018_TME"
# build_name="hg19"
# panel_name="TME_pan12"
# inAmp="/g/data3/yo4/Cancer-Epigenetics-Data/Nextseq500/wenjia/$project/config/TME_pan12_name_sorted.bed"
# add=200
# GREF="/g/data3/yo4/annotations/WGBS10X/${build_name}/${build_name}.fa"
# head $inAmp | column -t
# #chromA  startA     endA       nameA    primerA     Number.of.CpG.excluding.primer
# chr1     23882888   23882928   ID3      P1_ID3      1
# chr1     24469272   24469325   IL22RA1  P2_IL22RA1  1
# chr1     221054327  221054372  HLX      P2_HLX      4
# chr1     225837415  225837470  ENAH     P1_ENAH     1
# chr10    17274434   17274484   VIM      P2_VIM      1
# chr10    77156233   77156287   ZNF503   P2_ZNF503   3
# chr11    12026825   12026860   DKK3     P1_DKK3     1
# chr11    32455683   32455742   WT1      P2_WT1      2
# chr11    34622136   34622193   EHF      P2_EHF      2

### 1. Extend the regions by +-200 and check amplicon overlapping given a bed file of amplicon starts and ends.  
### If there is an overlapping, these overlapped amplicons are merged into one amplicon
bedtools merge -o collapse -delim "." -c 4 -i <( awk -v add="$add" '{OFS="\t"}NR>1{print $1, $2-add, $3+add, $5}' "$inAmp" | sort -k1,1 -k2,2n ) | \
sort -k1,1 -k2,2n > "$inAmp_extend_merged_bed"
# head "${inAmp/.bed/.extend.merged.bed}"| column -t
# chr1   23882688   23883128   P1_ID3
# chr1   24469072   24469525   P2_IL22RA1
# chr1   221054127  221054572  P2_HLX
# chr1   225837215  225837670  P1_ENAH
# chr10  17274234   17274684   P2_VIM
# chr10  77156033   77156487   P2_ZNF503
# chr11  12026625   12027060   P1_DKK3
# chr11  32455483   32455942   P2_WT1

### 2. Create a map of hg19 coordinate vs amplicon coordinate and a map of merged amplicons
bedtools intersect -a "$inAmp_extend_merged_bed" -b <(awk '{OFS="\t"}NR>1{print $1, $2, $3, $5}' "$inAmp" | sort -k1,1 -k2,2n) -wa -wb | \
awk '{OFS="\t"}{print $0, $6-$2, $7-$2, $7-$6}' | sort -k1,1 -k2,2n > "$inAmp_sorted_mapped"
# cat "${inAmp/_sorted.bed/_sorted_mapped.bed}" | column -t
# chr1   23882688   23883128   P1_ID3           chr1   23882888   23882928   P1_ID3           200  240  40
# chr1   24469072   24469525   P2_IL22RA1       chr1   24469272   24469325   P2_IL22RA1       200  253  53
# chr1   221054127  221054572  P2_HLX           chr1   221054327  221054372  P2_HLX           200  245  45
# chr1   225837215  225837670  P1_ENAH          chr1   225837415  225837470  P1_ENAH          200  255  55
# chr10  17274234   17274684   P2_VIM           chr10  17274434   17274484   P2_VIM           200  250  50
# chr10  77156033   77156487   P2_ZNF503        chr10  77156233   77156287   P2_ZNF503        200  254  54
# chr11  12026625   12027060   P1_DKK3          chr11  12026825   12026860   P1_DKK3          200  235  35
# chr11  32455483   32455942   P2_WT1           chr11  32455683   32455742   P2_WT1           200  259  59
# chr11  34621936   34622393   P2_EHF           chr11  34622136   34622193   P2_EHF           200  257  57
# chr11  67350776   67351238   P2_GSTP1         chr11  67350976   67351038   P2_GSTP1         200  262  62

### 3. Get sequences
build_name_panel="${build_name}_${panel_name}"
refs="$REF/${build_name_panel}.fa"
# GREF="/g/data3/yo4/annotations/WGBS10X/${build_name}/${build_name}.fa"
# REF="/g/data3/yo4/Cancer-Epigenetics-Data/Nextseq500/wenjia/$project/ref"
# bedtools getfasta -name -fi $ref -bed <(awk '{OFS="\t"}{print $1, $2, $3, $4"::"$1":"$2"-"$3}' $amp) -fo "$newRefp/hg19_TME_pan12.fa"
# bedtools getfasta -name -fi $GREF -bed <(awk '{OFS="\t"}{print $1, $2, $3, $4}' "$inAmp_extend_merged_bed") -fo "$refs"
mkdir -p "$REF/tmp"
rm $refs

for i in `awk '{OFS=";"}!/^#/{print $1,$2,$3,$5}' "$inAmp"`; do 
	chr=$(echo $i| cut -d";" -f1);
	s0=$(bc <<< "$(echo $i| cut -d";" -f2) - $add");
	s1=$(bc <<< "$(echo $i| cut -d";" -f2) - $add + 1"); 
	e=$(bc <<< "$(echo $i| cut -d";" -f3) + $add"); 
	n=$(echo $i| cut -d";" -f4); 
	fn="${n}::${chr}:${s0}-${e}"; 
	l="${chr}:${s1},${e}"; 
	wget -O "$REF/tmp/${fn}.fasta" "http://genome.ucsc.edu/cgi-bin/das/${build_name}/dna?segment=$l";
	sleep 3;
done

for i in `find "$REF/tmp" -name "*.fasta"`; do 
	echo $i; 
	echo -e ">$(basename ${i/.fasta/})" >> "$refs"; 
	awk '!/^</' $i| tr -d "\n" >> "$refs";
	echo -e "\n" >> "$refs" 
done
rm -r "$REF/tmp"

### 4. indexing with bismark
# bismark="/short/yo4/lpl913/software/Bismark/v0.22.1"
# module load bowtie2/2.2.5
# https://github.com/FelixKrueger/Bismark/blob/master/Bismark_User_Guide.md
"$bismark/bismark_genome_preparation" --verbose --bowtie2 "$REF"

### 5. Create CpG reference
# create_CpG_ref="/home/913/lpl913/allpipe/MethPanel/script/index_ref/create_CpG_ref.r"
echo "$build_name_panel"
`which R` -f "$create_CpG_ref" --args "$build_name_panel" "$add" "$REF/"

# make CpG mapped
frame="$REF/${build_name_panel}.CpG.bed"
bedtools intersect -a <(sort -k1,1 -k2,2n $frame) -b <(awk '{OFS="\t"}{print $4"::"$1":"$2"-"$3, $9, $10, $8"::"$5":"$6"-"$7}' "$inAmp_sorted_mapped"| sort -k1,1 -k2,2n) -loj | \
awk '{OFS="\t"}{i=1; if($8=="."){$8=$1; i=0}; print $1, $2, $3, i, $8}'| sort -k1,1 -k2,2n > "${frame/.bed/.mapped.bed}"

# check final
if [[ -f "$refs" ]]; then echo "done" > "$REF/indexing.txt"; fi


