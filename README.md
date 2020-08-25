# *MethPanel*

*MethPanel* is a computational pipeline in Linux operating system with an interactive graphical interface for rapid analysis of multiplex bisulphite PCR sequencing data. The tool covers a complete analysis workflow from genomic alignment to DNA methylation calling and supports an unlimited number of PCR amplicons and input samples. Moreover *MethPanel* offers important and unique features, such as a polymorphism score and a bisulphite PCR bias correction. *MethPanel* can be run in parallel by samples on either a personal computer or a high performance computer. The outputs are automatically forwarded to a shinyApp for convenient display, visualisation and sharing of data with collaborators and clinicians.

## *MethPanel* workflow
<img src="https://raw.githubusercontent.com/thinhong/MethPanel/master/figures/full_workflow_ver5.png" height="320">

### Installation
* *MethPanel* is built based on in-house bash/python/R script, Bpipe and a collection of software packages:
  * [bowtie2, 2.3.5.1]                   http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  * [Bismark, 0.22.3]                    https://www.bioinformatics.babraham.ac.uk/projects/bismark/
  * [bedtools, 2.28.0]                   https://github.com/arq5x/bedtools2
  * [samtools, 1.9]                      https://github.com/samtools/
  * [Bpipe, 0.9.9.9]                     https://github.com/ssadedin/bpipe
  * [fastqc, 0.11.9]                     http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  * [trim_galore, 0.6.5]                 https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
  * [UCSC-format-file-converter, 1.0]    http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
  * [BSgenome.Hsapiens.UCSC.version]     https://github.com/Przemol/BSgenome.Hsapiens.UCSC.version
  * [GenomicRanges]                      https://github.com/Bioconductor-mirror/GenomicRanges
  * [rtracklayer]                        https://github.com/Bioconductor-mirror/rtracklayer
  * [python, 3.7.4]                      https://www.python.org/
  * [python modules] configobj, argparse
  * [R, 3.6.1]                           https://www.r-project.org/
  * [R packages] ggplot2, rtracklayer, Biostrings, data.table, gplots, reshape2, shiny, shinyauthr, shinyjs, plotly
  * [multiqc]                            https://multiqc.info/
  * [java, jdk-13.33]                    https://java.com/en/download/
  * [shiny server] Visit our [**Wiki page**](https://github.com/thinhong/MethPanel/wiki/2.-shinyApp-client) for a detailed installation manual of shiny server and shinyapp
  
* Dependencies
  * All the above listed software packages, here are some commands to install python and R packages
    * [python version >= 3.7.4]
      ```
      python -m pip install configobj argparse
      ```
    * [R version >= 3.6.1] 
      ```
      source("https://bioconductor.org/biocLite.R")
      biocLite(c("BSgenome.Hsapiens.UCSC.hg19", "rtracklayer", "GenomicRanges")) 
      install.packages(c("shiny", "shinyauthr", "shinyjs", "plotly", "ggplot2", "gplots", "reshape2"))
      ```

* Installation
  ```
  git clone git@github.com:/thinhong/MethPanel.git
  ```
  * or
  ```
  git clone https://github.com/thinhong/MethPanel.git
  ```
  
* Inputs required
  * Single/paired fastq files in .gz format for each sample
  * Filled sample config file, [example](https://raw.githubusercontent.com/thinhong/MethPanel/master/config/sample.Example.pre.config)
  * Filled system config file, [example](https://raw.githubusercontent.com/thinhong/MethPanel/master/config/system.Example.pre.config)
  * DNA methylation marker panel file, [example](https://raw.githubusercontent.com/thinhong/MethPanel/master/config/amplicons.tsv)
  * Note: 
       * The fastq files are located in `/path/to/upstream/analysis/<user>/${project}/raw/Sample1/Sample1_R1.fastq.gz`,                                        `/path/to/upstream/analysis/<user>/${project}/raw/Sample1/Sample1_R2.fastq.gz` (if paired-end) for each sample.

### How to run *MethPanel*
```
project="<project_name>"
sample_config="/path/to/upstream/analysis/<user>/${project}/config/sample.${project}.pre.config"
system_config="/path/to/upstream/analysis/<user>/${project}/config/system.${project}.pre.config"
```
The fastq files are located in `/path/to/upstream/analysis/<user>/<project_name>/raw/`
```
python "/path/to/pipe/run_Bpipe.py" $sample_config $system_config
```
## *MethPanel* shinyApp
<img src="https://raw.githubusercontent.com/thinhong/MethPanel/master/figures/shiny_ver5.png" height="500">

For further details manual of *MethPanel* shinyApp, please visit our [**Wiki page**](https://github.com/thinhong/MethPanel/wiki). 
