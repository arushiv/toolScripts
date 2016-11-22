library(ggplot2)
## library(reshape2)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

colnames(d) <-  c("gene","snp","beta","tstat","pvalue")
mind <- aggregate(pvalue ~ gene, d, function(x) min(x))
mind <- subset(mind, pvalue<=0.05)

write.table(merge(mind, d), file=args[2], sep="\t", row.names=FALSE, quote=FALSE)

