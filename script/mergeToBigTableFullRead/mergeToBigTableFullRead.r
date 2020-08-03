# module load R/3.6.1
# get parameters
args = commandArgs(trailingOnly=TRUE)
inpath <- args[1]
outpath <- args[2]

### DNA meth for each CpG sites of amplicons
rmcols <- c("mean_ratio", "polymorphism", "Amplicon", "mean_cov", "Sample")
rmcols <- c("mean_ratio", "polymorphism", "mean_cov", "Sample")
fns <- list.files(inpath, pattern="*.polymorphism", full.name=T)
### DNA meth for each CpG sites
tb <- lapply(fns, function(x){
	dt <- read.table(x, header=T)
	tmeth <- setNames(data.frame(t(dt[,!(colnames(dt) %in% rmcols)])), paste0(dt[,"Sample"], ".ratio"))
	tcov <- setNames(data.frame(t(dt[,rep("mean_cov", dim(tmeth)[1])])), paste0(dt[,"Sample"], ".cov"))
	rownames(tcov) <- rownames(tmeth)
	t <- cbind(tcov, tmeth)
	t <- t[grep("Amplicon", rownames(t), invert=T),]
	t$Gpos <- as.integer(gsub("X", "", rownames(t)))
	t$Amplicon <- dt$Amplicon[1]
	return(t)
})

# btb <- do.call(rbind, tb)
df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all=T), tb)
df$Cpos <- df$Gpos - 2
colns <- sort(colnames(df))
colns <- c("Amplicon", "Cpos", "Gpos", colns[!(colns %in% c("Amplicon", "Cpos", "Gpos"))])
df <- df[colns]
df <- df[order(df$Amplicon, df$Gpos),]
# write raw matrix b to csv file
write.table(df, paste0(outpath, "/bigTable_fullRead.tsv"), quote=F, sep="\t", row.names=F, col.names=T)
# head(df)

### average DNA meth/polymorphism for whole of amplicon
keepcols <- c("mean_ratio", "mean_cov", "polymorphism")
tb <- lapply(fns, function(x){
	dt <- read.table(x, header=T)
	t <- setNames(data.frame(t(dt[, keepcols])), dt[,"Sample"])
	t$Amplicon <- dt$Amplicon[1]
	t$Type <- keepcols
	return(t)
})
# btb <- do.call(rbind, tb)
df <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all=T), tb)
colns <- sort(colnames(df))
colns <- c("Amplicon", "Type", colns[!(colns %in% c("Amplicon", "Type"))])
df <- df[colns]
df <- df[order(df$Amplicon, df$Type),]
# write raw matrix b to csv file
write.table(df, paste0(outpath, "/bigTable_fullRead_polymorphism.tsv"), quote=F, sep="\t", row.names=F, col.names=T)
# head(df)
