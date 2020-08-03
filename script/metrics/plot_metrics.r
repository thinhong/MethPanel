# module load R/3.6.1
# R
# library(knitr)
# library(kableExtra)
# library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
p <- args[1]


fn <- grep("bigTable_filtered.tsv.", list.files(p, pattern="*.tsv", full.names=T), value=T)
t <- read.table(fn, header=T)
coln <- colnames(t)
i <- which(coln=="Sample")
rn <- paste0(colnames(t), ".filtered")
rn[i] <- "Sample"
colnames(t) <- rn

fn <- grep("bigTable_fullRead.tsv.", list.files(p, pattern="*.tsv", full.names=T), value=T)
f <- read.table(fn, header=T)
coln <- colnames(f)
i <- which(coln=="Sample")
rn <- paste0(colnames(f), ".fullRead")
rn[i] <- "Sample"
colnames(f) <- rn

fn <- grep("metrics.tsv", list.files(p, pattern="*.tsv", full.names=T), value=T)
m <- read.table(fn, header=T)
coln <- colnames(m)
i <- which(coln=="Sample")
rn <- paste0(colnames(m), ".m")
rn[i] <- "Sample"
colnames(m) <- rn

tb <- Reduce(function(x, y){ merge(x, y, by='Sample', all=T) }, list(m, f, t))
cols <- c("Sample", "Number Read", "Mapped Read (%)", "Mapped converted (%)", "Mapped converted 40mapq (%)", "Methylation CpG (%)", "Methylation nonCpG (%)", "Depth of CpG (X)" , "Depth of CpG with full read (X)", "No. CpG dropout", "No. Amplicon dropout", "No. CpG dropout with full read", "No. Amplicon dropout with full read")
tb[,"Number Read"] <- tb[,"Number_Read.m"]
tb[,"Mapped Read (%)"] <- round(100*tb[,"Number_Mapped_Read.m"]/tb[,"Number_Read.m"], 2)
tb[,"Mapped converted (%)"] <- round(100*tb[,"Number_Mapped_Converted_Read.m"]/tb[,"Number_Read.m"], 2)
tb[,"Mapped converted 40mapq (%)"] <- round(100*tb[,"Number_Mapped_Converted_40mapq_Read.m"]/tb[,"Number_Read.m"], 2)
tb[,"Methylation CpG (%)"] <- gsub("%", "", tb[,"Methylation_CpG.m"])
tb[,"Methylation nonCpG (%)"] <- gsub("%", "", tb[,"Methylation_nonCpG.m"])
tb[,"Depth of CpG (X)"] <- round(tb[,"mean_cov.filtered"], 2)
tb[,"Depth of CpG with full read (X)"] <- round(tb[,"mean_cov.fullRead"], 2)
tb[,"No. CpG dropout"] <- tb[,"CpG_no_cov.filtered"]
tb[,"No. Amplicon dropout"] <- tb[,"Amp_no_cov.filtered"]
tb[,"No. CpG dropout with full read"] <- tb[,"CpG_no_cov.fullRead"]
tb[,"No. Amplicon dropout with full read"] <- tb[,"Amp_no_cov.fullRead"]

write.table(tb[cols], paste0(p, "metrics_table.tsv"), sep="\t", col.names=T, row.names=F)

# kable(tb[cols], align="c") %>%
# 	kable_styling("striped", full_width=T) %>%
# 	row_spec(0, align='c') %>%
# 	save_kable(file=paste0(p, "/metrics_table.html"), self_contained=T)
