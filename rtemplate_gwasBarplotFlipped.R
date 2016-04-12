library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stats)
library(qvalue)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

# trait   cell	motif	 enhancer_type    threshold       overlap expected_overlap        pval y_order
#levels(d$trait) <- c("Fast.Glucose(35)", "Fast.GluBMIadj(11)","Fast.Insulin(16)","Fast.InsulinBMIadj(15)","RandomGWASvariants(89)","Rheumatoid.Arthiritis(27)","T2D(95)")
d$trait <- factor(d$trait, levels=c("T2D","FGlu","T1D","Systemic.lupus.erythematosus"))
d$log2_enrichment <- log2(d$overlap/d$expected_overlap)
d$signedminuslog_pval <- -(log10(d$pval))*((d$log2_enrichment/abs(d$log2_enrichment)))


## pdf(paste("bubble.",args[2],sep=""), height=16, width=10)
## ## for (i in 1:length(vector)){

## ## maxE <- max(df$log2_enrichment, na.rm=T)
## ## minE <- min(df$log2_enrichment, na.rm=T)
## ggplot(d, aes(x=signedminuslog_pval, y=motif)) + geom_point(aes(fill=log2_enrichment, shape=cell), size=1.5) + scale_fill_gradientn(colours=c("white","red","darkred"), values=rescale(c(0,minE,maxE)), breaks=c(0,minE,maxE), labels=c(0,minE,maxE)) +  facet_wrap(~trait, nrow=1) + labs(title="Enrichment for GWAS variants for different traits in TF motif footprints in ATAC-seq", x="signed -log(p value)", y="motif") + theme(strip.text.x = element_text(size = 11), panel.background=element_rect(fill='white', colour='black'), axis.text.y=element_text(size=8)) + scale_shape_manual(values=c(21, 22, 23))
## dev.off()
pdf(args[2], height=12, width=9)
## ggplot(d, aes(x=motif, y=log2_enrichment)) + geom_bar(aes(fill=cell),stat="identity", position="dodge", width=0.8) + labs(title="Enrichment for GWAS variants for different traits in TF motif footprints in ATAC-seq", y="log2(fold enrichment)", y="motif") + facet_wrap(~trait,ncol=1, scales="free_x") + theme(strip.text.x = element_text(size = 8), panel.background=element_rect(fill='white', colour='black'), axis.text.y=element_text(size=6), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4))
ggplot(d, aes(x=motif, y=log2_enrichment)) + geom_bar(aes(fill=cell),stat="identity", position="dodge", width=0.8) + coord_flip() + labs(title="Enrichment for GWAS variants for different traits in TF motif footprints in ATAC-seq", y="log2(fold enrichment)", y="motif") + facet_wrap(~trait, nrow=1) + theme(strip.text.x = element_text(size = 8), panel.background=element_rect(fill='white', colour='black'), axis.text.y=element_text(size=10), axis.text.x=element_text(vjust=0.5, hjust=0.5, size=10), panel.grid.major.y=element_line(colour='black', size=0.08), panel.grid.major.x=element_blank())

## #density(d$pval)
## ## plot(hist(d$pval))
## ## }
dev.off()

