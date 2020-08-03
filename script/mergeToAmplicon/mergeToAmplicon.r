library(data.table)
args <- commandArgs(trailingOnly=TRUE)
inpath <- args[1]
amp <- args[2]
outpath <- args[3]

fns <- grep("/pattern/", list.files(inpath, recursive=T, pattern=paste0(amp, ".tsv"), full.names=T), value=T)
dt <- lapply(fns, function(x){
		Sample <- basename(dirname(dirname(x)));
		tb <- read.table(x, header=T, stringsAsFactors=F);
		tb$Sample <- Sample;
		tb$Amplicon <- amp;
		t <- t(tb[dim(tb)[1],]);
		t <- rbind(t, mean_cov=sum(tb$Frequency[-dim(tb)[1]]))
		colnames(t) <- Sample;
		return(t);
	})

tb <- t(Reduce(function(x, y) merge(x, y, all=TRUE), lapply(dt, function(y) data.table(y, keep.rownames=TRUE, key="rn"))))
rownames(tb) <- 1:dim(tb)[1]
colnames(tb) <- tb[1,]
tb <- tb[-1,]
tb <- as.data.frame(tb)
colX <- paste0("X", sort(as.integer(gsub("X", "", grep("^X", colnames(tb), value=T)))))
colL <- grep("^X", colnames(tb), value=T, invert=T)
tb <- tb[c(colX, colL)]
colnames(tb) <- gsub("Frequency", "mean_ratio", colnames(tb))
colnames(tb) <- gsub("Percentage", "polymorphism", colnames(tb))
colnames(tb) <- gsub("X", "", colnames(tb))
write.table(tb, paste0(outpath, "/", amp, ".polymorphism"), quote=F, sep="\t", row.names=F, col.names=T)
