# Clone MethPanel github
git clone git@github.com:/thinhong/MethPanel.git 

# Download software
mkdir -p MethPanel/software
cd MethPanel/software
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip
wget https://github.com/FelixKrueger/Bismark/archive/0.22.3.tar.gz
wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
wget https://github.com/ssadedin/bpipe/releases/download/0.9.9.9/bpipe-0.9.9.9.tar.gz
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.zip
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget https://github.com/ewels/MultiQC/archive/v1.9.zip

# install python modules with python version >= 3.7.4
pip install --user configobj argparse

# install R and Bioconductor package with R version >= 3.6.1
R -e '''
install.packages(c("ggplot2", "data.table", "gplots", "reshape2", "shiny", "shinyauthr", "shinyjs", "plotly"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "GenomicRanges", "rtracklayer"))
q()
'''
# please go to the directory (MethPanel/software) to install all the dependencies.
