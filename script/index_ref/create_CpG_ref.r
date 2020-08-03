### find CpG sites in R
# module load R/3.6.1
library(rtracklayer)
library(Biostrings)
options(scipen=999)

args = commandArgs(trailingOnly=TRUE)
build.name <- args[1]
add <- as.numeric(args[2])
outpath <- args[3]
fastap <- paste0(outpath, build.name, ".fa")

##############################################################################################################################
extractCytosinesFromFASTA <- function(file, contexts=c('CG','CHG','CHH'), anchor.C=NULL) {
    ### C anchors ###
    if (is.null(anchor.C)) {
        anchor.C <- regexpr('C', contexts)[1:length(contexts)]
        names(anchor.C) <- contexts
    }
    if (!all(names(anchor.C)==contexts)) {
        stop("names(anchor.C) must be equal to contexts")
    }
    ### Contexts ###
    if (is.null(names(contexts))) {
        context.names <- contexts
        dup.contexts <- duplicated(contexts) | duplicated(contexts, fromLast=TRUE)
        context.names[dup.contexts] <- paste0(contexts, '_', anchor.C)[dup.contexts]
    } else {
        context.names <- names(contexts)
    }

    ### Read file
    fasta <- Biostrings::readDNAStringSet(file)
    
    ### Chromosome lengths
    chromlengths <- width(fasta)
    names(chromlengths) <- names(fasta)
    
    ### Positions with ambiguous nucleotide codes
    notuse <- GRangesList()
    ambicodes <- c('R','Y','S','W','K','M','B','D','H','V','N')
    for (ambicode in ambicodes) {
        notuse.string <- Biostrings::vmatchPattern(ambicode, subject=fasta, fixed=TRUE)
        notuse[[ambicode]] <- as(notuse.string, 'GRanges')
    }
    notuse <- unlist(notuse, use.names=FALSE)
    
    ### Extract cytosine positions with context
    cytosines.forward <- GRangesList()
    for (icontext in 1:length(contexts)) {
        context <- contexts[icontext]
        context.name <- context.names[icontext]
        positions <- Biostrings::vmatchPattern(context, subject=fasta, fixed=FALSE)
        positions <- as(positions, 'GRanges')
        strand(positions) <- '+'
        positions$context <- factor(context.name, levels=unique(context.names))
        cytosines.forward[[icontext]] <- positions
    }
    cytosines.forward <- unlist(cytosines.forward, use.names=FALSE)
    seqlengths(cytosines.forward) <- chromlengths
    # Remove ambiguous positions
    ind <- findOverlaps(cytosines.forward, notuse)
    if (length(ind) > 0) {
        cytosines.forward <- cytosines.forward[-ind@from]
    }
    # Set width to 1
    end(cytosines.forward) <- start(cytosines.forward)
    
    ### Extract reverse context
    fasta <- Biostrings::reverseComplement(fasta)
    cytosines.reverse <- GRangesList()
    for (icontext in 1:length(contexts)) {
        context <- contexts[icontext]
        context.name <- context.names[icontext]
        positions <- Biostrings::vmatchPattern(context, subject=fasta, fixed=FALSE)
        positions <- as(positions, 'GRanges')
        strand(positions) <- '-'
        positions$context <- factor(context.name, levels=unique(context.names))
        cytosines.reverse[[icontext]] <- positions
    }
    cytosines.reverse <- unlist(cytosines.reverse, use.names=FALSE)
    seqlengths(cytosines.reverse) <- chromlengths
    # Correct coordinates from reverseComplement
    starts <- chromlengths[as.character(seqnames(cytosines.reverse))] - start(cytosines.reverse) + 1
    ends <- chromlengths[as.character(seqnames(cytosines.reverse))] - end(cytosines.reverse) + 1
    cytosines.reverse <- GRanges(seqnames=seqnames(cytosines.reverse), ranges=IRanges(start=ends, end=starts), strand='-', context=cytosines.reverse$context)
    # Remove ambiguous positions
    ind <- findOverlaps(cytosines.reverse, notuse)
    if (length(ind) > 0) {
        cytosines.reverse <- cytosines.reverse[-ind@from]
    }
    # Set width to 1
    start(cytosines.reverse) <- end(cytosines.reverse)
    
    ## Merge
    cytosines <- c(cytosines.forward, cytosines.reverse)
    
    # Shift positions by position of anchor C
    cind <- anchor.C[cytosines$context]
    strandint <- c('-'=-1,'+'=1,'*'=1)[as.character(strand(cytosines))]
    starts <- start(cytosines) + strandint * cind - strandint * 1
    ends <- end(cytosines) + strandint * cind - strandint * 1
    cytosines <- GRanges(seqnames=seqnames(cytosines), ranges=IRanges(start=starts, end=ends), strand=strand(cytosines), context=cytosines$context)
    seqlengths(cytosines) <- chromlengths
    
    ## Sort
    cytosines <- sort(cytosines, ignore.strand=TRUE)
    
    return(cytosines)
}
####################################################################################################################################################

### read fasta to sequence object
dna <- readDNAStringSet(fastap)

# CpG
CpG <- extractCytosinesFromFASTA(fastap, 'CG')
score(CpG) <- 0
for (id in grep("lambda", names(dna), value=T, invert=T)) {
	len <- width(dna[names(dna)==id]) - add
	ind <- (start(CpG[seqnames(CpG)==id])>=(add-1))&(end(CpG[seqnames(CpG)==id])<=(len+1))
	if(sum(ind)>0) { score(CpG[seqnames(CpG)==id][ind]) <- 1 }
}
sum(score(CpG)>0)
# [1] 316
# save to file as bed
export(CpG, paste0(outpath, build.name, ".CpG.strands.bed"), format="bed")

# CHH
CHH <- extractCytosinesFromFASTA(fastap, 'CHH')
score(CHH) <- 0
for (id in grep("lambda", names(dna), value=T, invert=T)) {
	len <- width(dna[names(dna)==id]) - add
	ind <- (start(CHH[seqnames(CHH)==id])>=add)&(end(CHH[seqnames(CHH)==id])<=len)
	if(sum(ind)>0) { score(CHH[seqnames(CHH)==id][ind]) <- 1 }
}
# save to file as bed
export(CHH, paste0(outpath, build.name, ".CHH.strands.bed"), format="bed")

# CHG
CHG <- extractCytosinesFromFASTA(fastap, 'CHG')
score(CHG) <- 0
for (id in grep("lambda", names(dna), value=T, invert=T)) {
	len <- width(dna[names(dna)==id]) - add
	ind <- (start(CHG[seqnames(CHG)==id])>=add)&(end(CHG[seqnames(CHG)==id])<=len)
	if(sum(ind)>0) { score(CHG[seqnames(CHG)==id][ind]) <- 1 }
}
# save to file as bed
export(CHG, paste0(outpath, build.name, ".CHG.strands.bed"), format="bed")

# CpG merged strands
CpGdf <- as.data.frame(CpG[which(strand(CpG)=="+")])
CpGdf$start <- CpGdf$start - 1
CpGdf$end <- CpGdf$end + 1
write.table(CpGdf[c("seqnames", "start", "end", "score")], paste0(outpath, build.name, ".CpG.bed"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# make a file which contains chromosomes size hg19.chrom.sizes.short
# get size of chromosomes
dt <- data.frame(chr=names(dna), len=width(dna))
# save
write.table(dt, paste0(outpath, build.name, ".chrom.sizes.short"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


