# module load R/3.6.1

library(ggplot2)
library(gplots)
library(reshape2)


colors <- colorRampPalette(c("lightblue", "green", "white", "yellow", "orange", "red"))(100)

args <- commandArgs(trailingOnly=TRUE)
inpath <- args[1]
suffix <- args[2]

amsp <- paste0(inpath, "bigTable", suffix, ".tsv")
metricp <- paste0(inpath, "metrics/bigTable", suffix, ".tsv")
system(paste0("mkdir -p ", dirname(metricp)))

### All
ams <- read.table(amsp, header=TRUE, sep="\t", comment.char = "\\", stringsAsFactors=F, check.names=FALSE)
cov <- grep(".cov", colnames(ams), value=T)
t <- ams[, cov]
t[t=="."] <- NA
t <- as.data.frame(sapply(t, as.integer))
CpG_no_cov <- paste0(colSums(apply(t, 2, function(x) is.na(x))), "/", dim(t)[1])

cov_ave <- aggregate(t, list(ams[,"Amplicon"]), function(x){round(mean(x, na.rm=T), 0)})
colname <- gsub("\\.cov", "", grep("\\.cov", colnames(cov_ave), value=T))
rowname <- cov_ave$Group.1
Amp_no_cov <- paste0(colSums(apply(cov_ave[,-1], 2, function(x){is.na(x)})), "/", dim(cov_ave)[1])
# plot average coverage of each amplicon
pdf(paste0(amsp, ".coverage.heatmap.pdf"))
heatmap.2(as.matrix(cov_ave[cov]), col=colors, dendrogram="none", cexRow=0.5, cexCol=0.5, labCol=colname, labRow=rowname,
	symkey=FALSE, density.info="none", trace="none", Colv=F, Rowv=F, cellnote=cov_ave[cov], notecex=0.35, 
	notecol="black", margins=c(12,5))
heatmap.2(as.matrix(cov_ave[cov]), col=colors, dendrogram="none", cexRow=0.5, cexCol=0.5, labCol=colname, labRow=rowname, scale="row",
	symkey=FALSE, density.info="none", trace="none", Colv=F, Rowv=F, cellnote=cov_ave[cov], notecex=0.35, 
	notecol="black", margins=c(12,5))
heatmap.2(as.matrix(cov_ave[cov]), col=colors, dendrogram="none", cexRow=0.5, cexCol=0.5, labCol=colname, labRow=rowname, scale="column",
	symkey=FALSE, density.info="none", trace="none", Colv=F, Rowv=F, cellnote=cov_ave[cov], notecex=0.35, 
	notecol="black", margins=c(12,5))
dev.off()

# plot boxplot
a <- melt(cov_ave)
a$Amp <- sapply(a$Group.1, function(x){ strsplit(x, "::")[[1]][1] })
ggplot(a, aes(y=value, x=Amp, color=Amp)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25) + xlab("Amplicons") + ylab("Coverage") +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = "none") +
ggsave(paste0(amsp, ".coverage.boxplot1.pdf"))
ggplot(a, aes(y=value, x=Amp, color=Amp)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Amplicons") + ylab("Coverage") +
theme(legend.position = "none") + ggsave(paste0(amsp, ".coverage.boxplot2.pdf"))

# make a csv file
write.table(cov_ave, paste0(amsp, ".coverage.tsv"), sep="\t", row.names=F, col.names=T)
mean_cov <- colMeans(cov_ave[,-1], na.rm=T)

# plot coverage of all amplicons for each sample
p <- list()
for(i in cov){
	print(i)
	a <- cov_ave[, c("Group.1", i)]
	a$Amp <- sapply(a$Group.1, function(x){ strsplit(x, "::")[[1]][1] })
	colnames(a) <- c("X", "cov", "Amp")
	p[[i]] <- ggplot(a, aes(y=cov, x=Amp)) + geom_point() + theme(axis.text.x=element_text(angle=90, hjust=1)) + theme(legend.position = "none") + ggtitle(i) + xlab("Amplicons") + ylab("Coverage")
	# ggsave(paste0(inpath, "bigTable.sorted.tsv.coverage.lineplot.", i, ".pdf"))
}
# ggsave(paste0(inpath, "bigTable.sorted.tsv.coverage.lineplot.eachsample.pdf"))
pdf(paste0(amsp, ".coverage.point.eachsample.pdf"), onefile = TRUE)
for(i in cov){
		print(p[[i]])
}
dev.off()

### DNA methylation
ratio <- grep(".ratio", colnames(ams), value=T)
t <- ams[, ratio]
t[t=="."] <- NA
t <- as.data.frame(sapply(t, as.integer))
ratio_ave <- aggregate(t, list(ams[,"Amplicon"]), function(x){round(mean(x),1)})

colname <- gsub(".ratio", "", ratio)
rowname <- cov_ave$Group.1
colors <- colorRampPalette(c("darkblue", "lightblue", "green", "yellow", "orange", "red"))(100)
pdf(paste0(amsp, ".methylation.heatmap.pdf"))
heatmap.2(as.matrix(ratio_ave[ratio]), col=colors, dendrogram="none", cexRow=0.5, cexCol=0.5, labCol=colname, labRow=rowname,
	symkey=FALSE, density.info="none", trace="none", Colv=F, Rowv=F, cellnote=ratio_ave[ratio], notecex=0.45, 
	notecol="black", margins=c(12,5))
heatmap.2(as.matrix(ratio_ave[ratio]), col=colors, dendrogram="none", cexRow=0.5, cexCol=0.5, labCol=colname, labRow=rowname, scale="row",
	symkey=FALSE, density.info="none", trace="none", Colv=F, Rowv=F, cellnote=ratio_ave[ratio], notecex=0.45, 
	notecol="black", margins=c(12,5))
heatmap.2(as.matrix(ratio_ave[ratio]), col=colors, dendrogram="none", cexRow=0.5, cexCol=0.5, labCol=colname, labRow=rowname, scale="column",
	symkey=FALSE, density.info="none", trace="none", Colv=F, Rowv=F, cellnote=ratio_ave[ratio], notecex=0.45, 
	notecol="black", margins=c(12,5))
dev.off()

# make a csv file
write.table(ratio_ave, paste0(amsp, ".methylation.csv"), sep="\t", row.names=F, col.names=T)
mean_ratio <- colMeans(ratio_ave[,-1], na.rm=T)
metric <- as.data.frame(t(rbind(mean_cov, mean_ratio, CpG_no_cov, Amp_no_cov)))
metric$Sample <- gsub(".cov", "", rownames(metric))
write.table(metric, paste0(metricp, ".coverage.methylation.tsv"), sep="\t", row.names=F, col.names=T)

# plot boxplot
a <- melt(ratio_ave)
a$Amp <- sapply(a$Group.1, function(x){ strsplit(x, "::")[[1]][1] })
ggplot(a, aes(y=value, x=Amp, color=Amp)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25) + xlab("Amplicons") + ylab("Methylation (%)") +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = "none") +
ggsave(paste0(amsp, ".methylation.boxplot1.pdf"))
ggplot(a, aes(y=value, x=Amp, color=Amp)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Amplicons") + ylab("Methylation (%)") +
theme(legend.position = "none") + ggsave(paste0(amsp, ".methylation.boxplot2.pdf"))

# plot methylation of all amplicons for each sample
p <- list()
for(i in ratio){
	print(i)
	a <- ratio_ave[, c("Group.1", i)]
	a$Amp <- sapply(a$Group.1, function(x){ strsplit(x, "::")[[1]][1] })
	colnames(a) <- c("X", "ratio", "Amp")
	p[[i]] <- ggplot(a, aes(y=ratio, x=Amp)) + geom_point() + theme(axis.text.x=element_text(angle=90, hjust=1)) + theme(legend.position="none") + ggtitle(i) + xlab("Amplicons") + ylab("Methylation (%)")
	# ggsave(paste0(inpath, "bigTable.sorted.tsv.coverage.lineplot.", i, ".pdf"))
}
# ggsave(paste0(inpath, "bigTable.sorted.tsv.coverage.lineplot.eachsample.pdf"))
pdf(paste0(amsp, ".methylation.point.eachsample.pdf"), onefile=TRUE)
for(i in ratio){
		print(p[[i]])
}
dev.off()

