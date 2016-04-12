library(ggplot2)
library(reshape2)

d <- read.table("", header=T)

pdf(".pdf", height=10, width=8)
qplot(data=d, x = id, y = values, geom="bar", stat = "identity") + facet_wrap("variable", scales="free") + theme(text = element_text(size=14), axis.text.x = element_text(angle=90, vjust=0.5, hjust=0.5))
dev.off()

vector <- levels(unique(df$enhancer_type))

pdf(args[2], height=8, width=20)
for (i in 1:length(vector)){
maxE <- max(subset(df,enhancer_type %in% c(vector[i]))$log2_enrichment, na.rm=T)
minE <- min(subset(df,enhancer_type %in% c(vector[i]))$log2_enrichment, na.rm=T)
print(ggplot(subset(df,enhancer_type %in% c(vector[i])), aes(x=signedminuslog_pval, y=cell), cex.lab=1) + geom_point(aes(fill=log2_enrichment, shape=threshold), size=2) + scale_fill_gradientn(colours=c("blue","white","red"), values=rescale(c(minE,0,maxE))) +  facet_wrap(~trait, nrow=1) + theme_bw() + labs(title=paste("Enrichment for GWAS variants for different traits in", vector[i], sep=" "), x="signed -log(p value)", y="cell type") + theme(strip.text.x = element_text(size = 11)) + geom_vline(xintercept=c(min(na.omit(subset(df,enhancer_type %in% c(vector[i]))$signedminuslog_pval)),-log10(0.05/1324)), linetype="dashed") + geom_rect(xmin=min(na.omit(subset(df,enhancer_type %in% c(vector[i]))$signedminuslog_pval)), xmax=-log10(0.05/1324), ymin = -Inf , ymax = Inf, alpha=0.002) + scale_shape_manual(values=c(21, 22, 23)))
}
dev.off()
