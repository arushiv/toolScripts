library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stats)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

# trait   cell		 enhancer_type    threshold       overlap expected_overlap        pval y_order
#levels(d$trait) <- c("Fast.Glucose(35)", "Fast.GluBMIadj(11)","Fast.Insulin(16)","Fast.InsulinBMIadj(15)","RandomGWASvariants(89)","Rheumatoid.Arthiritis(27)","T2D(95)")

d$log2_enrichment <- log2(d$overlap/d$expected_overlap)
d$signedminuslog_pval <- -(log10(d$pval))*((d$log2_enrichment/abs(d$log2_enrichment)))
df <- transform(d, cell=reorder(cell, -y_order))

#maxE <- max(df$log2_enrichment, na.rm=T)
#minE <- min(df$log2_enrichment, na.rm=T)


vector <- levels(unique(df$enhancer_type))

pdf(args[2], height=8, width=20)
for (i in 1:length(vector)){
maxE <- max(subset(df,enhancer_type %in% c(vector[i]))$log2_enrichment, na.rm=T)
minE <- min(subset(df,enhancer_type %in% c(vector[i]))$log2_enrichment, na.rm=T)
print(ggplot(subset(df,enhancer_type %in% c(vector[i])), aes(x=signedminuslog_pval, y=cell), cex.lab=1) + geom_point(aes(fill=log2_enrichment, shape=threshold), size=2) + scale_fill_gradientn(colours=c("blue","white","red"), values=rescale(c(minE,0,maxE))) +  facet_wrap(~trait, nrow=1) + theme_bw() + labs(title=paste("Enrichment for GWAS variants for different traits in", vector[i], sep=" "), x="signed -log(p value)", y="cell type") + theme(strip.text.x = element_text(size = 11)) + geom_vline(xintercept=c(min(na.omit(subset(df,enhancer_type %in% c(vector[i]))$signedminuslog_pval)),-log10(0.05/1324)), linetype="dashed") + geom_rect(xmin=min(na.omit(subset(df,enhancer_type %in% c(vector[i]))$signedminuslog_pval)), xmax=-log10(0.05/1324), ymin = -Inf , ymax = Inf, alpha=0.002) + scale_shape_manual(values=c(21, 22, 23)))
}
dev.off()



#pdf(args[2], height=8, width=20)
#ggplot(df, aes(x=signedminuslog_pval, y=cell), cex.lab=1) + geom_point(aes(fill=log2_enrichment, shape=threshold), size=2) + scale_fill_gradientn(colours=c("blue","white","red"), values=rescale(c(minE,0,maxE))) +  facet_wrap(~trait, nrow=1) + theme_bw() + labs(title="Enrichment for GWAS variants for different traits in enhancer state segmentations", x="signed -log(p value)", y="cell type") + theme(strip.text.x = element_text(size = 11)) + geom_vline(xintercept=c(min(na.omit(d$signedminuslog_pval)),-log10(0.05/1324)), linetype="dashed") + geom_rect(xmin=min(na.omit(d$signedminuslog_pval)), xmax=-log10(0.05/1324), ymin = -Inf , ymax = Inf, alpha=0.002) + scale_shape_manual(values=c(21, 22, 23))
#dev.off()


# pdf(args[2], height=8, width=20)
# ggplot(df, aes(x=log2_enrichment, y=cell), cex.lab=1) + geom_point(aes(colour=signed_minuslog_pval, shape=threshold)) + scale_colour_gradientn(colours=c("blue","white","red"), values=rescale(c(minP,0,maxP))) +  facet_wrap(~trait, nrow=1, scales="free_x") + theme_bw() + ggtitle("signed -log(p values) and log2(enrichment)")
# dev.off()

# pdf(args[2], height=8, width=20)
# ggplot(df, aes(x=signed_minuslog_pval, y=cell), cex.lab=1) + geom_point(aes(colour=log2_enrichment, shape=threshold)) + scale_colour_gradientn(colours=c("blue","white","red"), values=rescale(c(minE,0,maxE))) +  facet_wrap(~trait, nrow=1, scales="free_x") + theme_bw() + ggtitle("signed -log(p values) and log2(enrichment)")
# dev.off()


# pdf(args[3], height=8, width=20)
# ggplot(df, aes(x=filt_signed_minuslog_corrected_pval, y=cell), cex.lab=1) + geom_point(aes(colour=filt_log2_enrich, shape=threshold)) + scale_colour_gradientn(colours=c("blue","white","red"), values=rescale(c(minE1,0,maxE1))) +  facet_wrap(~trait, nrow=1) + theme_bw() + ggtitle("signed -log(p values) and log2(enrichment) filtered after Bonferroni correction")
# dev.off()



# pdf(args[4], height=8, width=20)
# ggplot(df, aes(x=signed_minuslog_corrected_pval, y=cell), cex.lab=1) + geom_point(aes(colour=log2_enrichment, shape=threshold)) + scale_colour_gradientn(colours=c("blue","white","red"), values=rescale(c(minE,0,maxE))) +  facet_wrap(~trait, nrow=1, scales="free_x") + theme_bw() + ggtitle("signed -log(bonferroni corrected p values) and log2(enrichment)")
# dev.off()
