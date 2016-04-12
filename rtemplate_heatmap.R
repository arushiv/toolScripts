library(ggplot2)
library(reshape2)
library(scales)
args <- commandArgs(TRUE)
df <- read.table(args[1], header=F)

df[,7] <- df[,6]*854
df[,8] <- -log(df[,6])
df[,9] <- -log(df[,7])
df[,10] <- log(df[,4] / df[,5])
df[,11] <- log2(df[,4] / df[,5])

colnames(df) <- c("trait","enhancer_type","cell","number","score","pval","corrected_pval","minus_log_pval","minus_log_corrected_pval","enrichment_values_log10", "enrichment_values_log2")

df

maxE <- max(df$enrichment_values_log2, na.rm=T)
minE <- min(df$enrichment_values_log2, na.rm=T)

maxP <- max(df$minus_log_corrected_pval, na.rm=T)
minP <- min(df$minus_log_corrected_pval, na.rm=T)
maxE
minE
maxP
minP

pdf(args[2], height=5, width=9)
ggplot(df, aes(x = cell, y = trait)) + geom_tile(aes(fill=minus_log_corrected_pval)) + scale_fill_gradientn(colours=c("grey","white", "pink", "red", "darkred"), values=rescale(c(0,1.3,7,15,maxP)), name="-log(Bonferroni corrected p value)") + facet_wrap("enhancer_type", scales="free") + theme(text = element_text(size=5), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), panel.background = element_rect(fill="white")) + geom_text(aes(label=number), size=1) + scale_colour_hue(l=10, c=250)
dev.off()


pdf(args[3], height=6, width=6)
	for (i in 1:length(vector)){
	maxE <- max(subset(d1,footprint %in% c(vector[i]))$enrich, na.rm=T)
	minE <- min(subset(d1,footprint %in% c(vector[i]))$enrich, na.rm=T)
	color_scale <- round(c(minE, 0, 0.33*maxE, .66*maxE, maxE), digits=2)
	print(ggplot(subset(d1, footprint %in% c(vector[i])), aes(x = chromatin_state, y = state_tissue)) + facet_grid(type~region) + geom_tile(aes(fill=enrich)) + scale_fill_gradientn(colours=c("lightblue","white","pink", "red", "darkred"), values=rescale(color_scale), breaks=color_scale, labels=color_scale, name="enrichment=log2(real_stat/mean_null)") + theme(text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y = element_text(size=4), panel.background = element_rect(fill="white")) + labs(title=paste("Enrichments for chromatin states overlaps with ATAC-seq",vector[i], sep=" "), x="Chromatin States", y="Cell/Tissue type") + scale_x_discrete(limits=x_order))
	}
	dev.off()