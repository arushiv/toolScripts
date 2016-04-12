library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stats)
library(qvalue)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

## eqtl    eqtl_fdr        quantile        bin_number      tissue  motif   overlap expected_overlap        pval
d <- subset(subset(d, quantile=="binQuintile"), overlap>=as.numeric(args[2]))

d$log2_enrichment <- log2(d$overlap/d$expected_overlap)
d$signedminuslog_pval <- -(log10(d$pval))*((d$log2_enrichment/abs(d$log2_enrichment)))
d$minuslog_pval <- -(log10(d$pval))
tomatch <- c("RFX","MA0509")
## Volcano plot
pdf(args[3], height=8, width=10)
ggplot(d, aes(y=minuslog_pval, x=log2_enrichment)) + geom_point(aes(fill=tissue), shape=21, size=2.5, alpha=0.8) + labs(title="log2(Fold enrichment) vs -log10(p value) for islet eQTLs in TF motif footprints in ATAC-seq", x="log2(fold enrichment)", y="-log10(p value)") + facet_grid(eqtl_fdr~bin_number) + theme(strip.text.x = element_text(size = 8), panel.background=element_rect(fill='white', colour='black'), axis.text.y=element_text(size=10), axis.text.x=element_text(vjust=0.5, hjust=0.5, size=10)) + geom_hline(aes(yintercept=4.75), colour='black', size=0.05) + geom_vline(aes(xintercept=0), colour='black', size=0.05) + geom_text(data=subset(d[grep(paste(tomatch,collapse="|"), d$motif),],minuslog_pval>=4.75), aes(label=motif), size=1, colour='black', hjust=0, vjust=0)
                                                                                                                                                                                                                                                                                             #panel.grid.major.y=element_line(colour='black', size=0.08), panel.grid.major.x=element_blank())
dev.off()

