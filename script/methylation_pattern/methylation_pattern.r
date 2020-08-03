# qsub -I -l mem=64G -l ncpus=8 -q expressbw -l walltime=24:00:00
# module load R/3.6.1
# R

library(data.table)
library(gplots)
colors <- colorRampPalette(c("lightblue", "green", "white", "yellow", "orange", "red"))(100)

args <- commandArgs(trailingOnly=TRUE)
inpath <- args[1]
outpath1 <- args[2]
outpath2 <- args[3]
outpath3 <- args[4]
sample <- args[5]
ampp <- args[6]
CpGp <- args[7]
CpG <- read.table(CpGp)

### read data
print(paste0("Read ", inpath))
dt <- fread(inpath, header=F)
amptb <- read.table(ampp, header=F)

# for each amplicon, do

for(j in 1:dim(amptb)[1]){
	# j <- 10 28 (2); 52 (1)
	amp <- paste0(amptb$V4[j], "::", amptb$V1[j], ":", amptb$V2[j], "-", amptb$V3[j])
	print(paste0("Amplicon == ", amp))
	# CpG in full amplicon
	CpGamp <- CpG[(CpG$V1==amp)&(CpG$V4==1), "V2"]
	CpGamp_strand <- unique(c(CpGamp, CpGamp + 1, CpGamp + 2))
	al <- as.data.frame(dt[dt$V3==amp, ])
	al <- al[al$V4 %in% CpGamp_strand,]
	# head(al)
	# a$V4 <- as.character(a$V4)
	if(dim(al)[1]>0){
		ampsub <- paste0(amptb$V8[j], "::", amptb$V5[j], ":", amptb$V6[j], "-", amptb$V7[j])
		aw <- dcast(al, V1 ~ V4, value.var="V5")

		# # find all the CpG sites that reads mapped to
		raw <- sort(as.integer(grep("V1", colnames(aw), value=T, invert=T)))
		# merge + and - strand
		ind <- which(diff(raw)==1)
		while(length(ind)>0){
			C <- as.character(raw[ind[1]])
			G <- as.character(raw[ind[1]+1])
			aw[,C] <- sub("NA", "", paste0(aw[, C], aw[, G]))
			aw <- aw[, !(names(aw) %in% G)]
			raw <- sort(as.integer(grep("V1", colnames(aw), value=T, invert=T)))
			ind <- which(diff(raw)==1)
		}
		raw <- sort(as.integer(grep("V1", colnames(aw), value=T, invert=T)))
		# replace NA by X
		aw[is.na(aw)] <- "X"
		aw[aw=="NA"] <- "X"

		startc <- as.integer(amptb[j, "V9"])
		# real ending of amplicon
		endc <- as.integer(amptb[j, "V10"])
		# amp_width
		amp_width <- endc - startc
		## Only keep CpG sites which are inside the amplicon width
		raw <- raw[(raw>=startc)&(raw<=endc)]
		#### find the frequecy of the read patterns
		if(length(raw)>=1){
			Xper <- vector("list", length(raw))
			names(Xper) <- raw
			for(i in raw){ x <- table(aw[,as.character(i)])['X']/sum(table(aw[,as.character(i)])); if(is.na(x)){x=0}; Xper[as.character(i)] <- x }
			raw_noX <- as.integer(names(which(Xper < 0.90)))
			if(length(raw_noX)>1){		
				g <- data.frame(table(apply(aw[,as.character(raw_noX)], 1, function(x) paste0(x, collapse=""))))
			}else if(length(raw_noX)==1){
				g <- data.frame(table(aw[, as.character(raw_noX)]))
			}else{
				print("No CpG in this amplicon after remove X")
				next
			}
		}else{
			print("No CpG in this amplicon before remove X")
			next
		}
		# g
		#  Var1 Freq
		# 1   ZZ   29
		# 2   Zz   28
		# 3   zZ    8
		# 4   zz   29
		g$Var1 <- gsub("z", "T", g$Var1)
		g$Var1 <- gsub("Z", "C", g$Var1)
		g <- g[order(g$Freq), ]
		# g
		#  Var1 Freq
		# 3   TC    8
		# 2   CT   28
		# 1   CC   29
		# 4   TT   29
		# write raw matrix b to csv file
		write.table(aw[,as.character(raw_noX)], paste0(outpath1, "/", ampsub, ".tsv"), quote=F, sep="\t", row.names=F, col.names=T)
		# make frequency table
		freq <- as.data.frame(t(apply(g, 1, function(x){c(strsplit(x["Var1"],"")[[1]], x["Freq"])})), stringsAsFactors=F)
		colnames(freq) <- c(raw_noX - startc, "Frequency")
		# freq
		# 2 41 Frequency
		# 3 T  C         8
		# 2 C  T        28
		# 1 C  C        29
		# 4 T  T        29
		write.table(freq, paste0(outpath2, "/", ampsub, ".tsv"), quote=F, sep="\t", row.names=F, col.names=T)
		### plotting
		freq <- freq[!apply(freq, 1, function(x){any(x=="X")}),]
		if(is.null(dim(freq)[1])||(dim(freq)[1]==0)){
			print("     + All reads contain Unknown!")
		}else if((dim(freq)[1]>=1)&&(dim(freq)[2]>=2)){
			freq$Frequency <- as.integer(freq$Frequency)
			colname <- as.character(raw_noX - startc)
			freq <- aggregate(freq['Frequency'], by=freq[colname], sum)
			freq <- freq[order(freq$Frequency),]
			freq["Percentage"] <- round(100*freq["Frequency"]/sum(freq["Frequency"]), 2)
			if(dim(freq)[1]==1){
				print("     + Only one pattern!")
				p <- 0
				freq <- data.frame(rbind(freq, freq))
				tb <- freq
				tb[tb=='C'] <- 100
				tb[tb=='T'] <- 0
				tb[tb=='X'] <- 50
				for(i in colname) tb[,paste0("X",i)] <- as.integer(tb[,paste0("X",i)])
				tb[dim(tb)[1]+1,] <- apply(tb, 2, function(x){ round(100*sum((x/100)*tb$Frequency)/sum(tb$Frequency),2) })
				tb[dim(tb)[1], "Frequency"] <- round(sum(tb[dim(tb)[1], paste0("X", colname)])/length(colname), 2)
				tb[dim(tb)[1], "Percentage"] <- 100*p
				freq[dim(tb)[1],] <- tb[dim(tb)[1],]
				write.table(freq, paste0(outpath3, "/", ampsub, ".tsv"), quote=F, row.names=F, sep="\t")
				pdf(paste0(outpath3, "/", ampsub, ".pdf"), width=10, height=10)
				heatmap.2(as.matrix(tb), col=colors, dendrogram="none", cex.main=1, scale="column", cexCol=1, main=paste0(ampsub, "\nLength=", amp_width, ", Sample = ", sample), symkey=T, 
				density.info="none", trace="none", Colv=F, Rowv=F, labCol=c(colname, "Frequency", "Percentage"), labRow="", cellnote=freq, notecex=1, notecol="black", margins=c(12,5), keysize=0.8, key.title="", key.xlab="")
				dev.off()
			}else{
				print("     + More than one pattern!")
				p <- round(1 - sum((freq$Percentage/100)^2), 4)
				tb <- freq
				tb[tb=='C'] <- 100
				tb[tb=='T'] <- 0
				tb[tb=='X'] <- 50
				for(i in colname) tb[,i] <- as.integer(tb[,i])
				tb[dim(tb)[1]+1,] <- apply(tb, 2, function(x){ round(100*sum((x/100)*tb$Frequency)/sum(tb$Frequency),2) })
				tb[dim(tb)[1], "Frequency"] <- round(sum(tb[dim(tb)[1], colname])/length(colname), 2)
				tb[dim(tb)[1], "Percentage"] <- 100*p
				freq[dim(tb)[1],] <- tb[dim(tb)[1],]
				write.table(freq, paste0(outpath3, "/", ampsub, ".tsv"), quote=F, row.names=F, sep="\t")
				pdf(paste0(outpath3, "/", ampsub, ".pdf"), width=10, height=10)
				heatmap.2(as.matrix(tb), col=colors, dendrogram="none", cex.main=1, scale="column", cexCol=1, main=paste0(ampsub, "\nLength=", amp_width, ", Sample = ", sample), symkey=T, 
				density.info="none", trace="none", Colv=F, Rowv=F, labCol=c(colname, "Frequency", "Percentage"), labRow="", cellnote=freq, notecex=1, notecol="black", margins=c(12,5), keysize=0.8, key.title="", key.xlab="")
				dev.off()
			}
		}
	}
}

