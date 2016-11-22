library(ggplot2)
library(reshape2)
args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

vector <- unique(subset(d,eQTLnum=='2')$GeneID)
d2 <- d[FALSE,]

d2 <- do.call(rbind, lapply(vector, function(i) {
    subset(d, (GeneID==i & eQTLnum=='2') | (GeneID==i & eQTLnum=='1'))
}))


d2 <- d2[-grep("HLA-", d2$GeneName),]

d2prim <- subset(d2, eQTLnum=='1')
d2sec <- subset(d2, eQTLnum=='2')
## d2[-grep(paste(tomatch, collapse="|"), d2$GeneID),]

snpinfo_prim <-  unique(paste("chr", d2prim$SNPchr, ":", d2prim$EndSNP, sep=""))
snpinfo_sec <-  unique(paste("chr", d2sec$SNPchr, ":", d2sec$EndSNP, sep=""))


write.table(d2, file=paste(args[2],"primSecGenes","dat",sep="."), sep='\t', quote=FALSE, row.names=FALSE)
write.table(snpinfo_prim, file=paste(args[2],"primary","fullset","txt",sep="."), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(snpinfo_sec, file=paste(args[2],"secondary","fullset","txt",sep="."), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
