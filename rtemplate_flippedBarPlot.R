library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stats)
library(qvalue)
args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)

# eqtl    eqtl_set        quantile        bin_number      feature cell    subfeature      overlap expected_overlap        pval    tissue

d$log2_enrichment <- log2(d$overlap/d$expected_overlap)
d$signedminuslog_pval <- -(log10(d$pval))*((d$log2_enrichment/abs(d$log2_enrichment)))
d$minuslog_pval <- -(log10(d$pval))
denhancerCluster <- subset(d, feature=="enhancerclusters")
denhancerCluster$cell <- factor(denhancerCluster$cell, levels=rev(c("cluster_5","cluster_27","cluster_4","cluster_41","cluster_43","cluster_56","cluster_2","cluster_42","cluster_13","cluster_28","cluster_31","cluster_29","cluster_32","cluster_36","cluster_22","cluster_54","cluster_7","cluster_46","cluster_40","cluster_47","cluster_12","cluster_19","cluster_58","cluster_16","cluster_60","cluster_57","cluster_51","cluster_59","cluster_11","cluster_9","cluster_10","cluster_30","cluster_37","cluster_38","cluster_49","cluster_24","cluster_39","cluster_53","cluster_45","cluster_21","cluster_48","cluster_1","cluster_35","cluster_26","cluster_52","cluster_20","cluster_3","cluster_25","cluster_44","cluster_55","cluster_18","cluster_6","cluster_33","cluster_15","cluster_17","cluster_8","cluster_34","cluster_50","cluster_23","cluster_14")))
## denhancerCluster$bin_number <- factor(denhancerCluster$bin_number, levels=c("NA","1","2","3","4","5"))
## qobj <- qvalue(d3$pval, fdr.level=0.05)
## d3$qsignificant <- qobj$significant
## d3$qvalue <- qobj$qvalues
## d3[,c("pval","qvalue")]

# Plot enrichment heatmap for fulleQTLset
## denhancerCluster_fullset <- subset(denhancerCluster, quantile=="fulleQTLset")
denhancerCluster$quantile <- factor(denhancerCluster$quantile, levels=c("fulleQTLset","binQuintile"))
pdf(paste(args[2],"enhancerCluster.pdf", sep="."), height=8, width=6)
lapply(sort(unique(denhancerCluster$eqtl_set)), function(i){
    dplot <- subset(denhancerCluster, eqtl_set==i)
    ## dplotsignificant <- subset(dplot, minuslog <- pval>=2.79)
    colVector <- c(min(dplot$log2_enrichment), 0, 0.35*(max(dplot$log2_enrichment)), 0.7*(max(dplot$log2_enrichment)), max(dplot$log2_enrichment))
    print(c(round(colVector,2)))
    ggplot(dplot, aes(y = signedminuslog_pval, x = cell)) + geom_bar(aes(fill=log2_enrichment), stat="identity", colour='black') + coord_flip()+ facet_wrap(quantile~bin_number, nrow=1) + theme(plot.title=element_text(size=5), axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=6), panel.background = element_rect(fill="white"), legend.position="bottom", legend.justification="center", legend.title=element_text(size=5),legend.text=element_text(size=4), legend.key.size=unit(6,"mm") , panel.grid.major.y=element_line(size=0.3, colour='grey'), panel.grid.major.x=element_blank()) + labs(title=paste("eQTL enrichment in enhancer state clusters",i,sep=" "), x="cluster", y="-log10(p value)") + scale_fill_gradientn(name="log2(fold enrichment)",colours=c("lightblue","white","pink","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_hline(yintercept=c(min(na.omit(dplot$signedminuslog_pval)),-log10(0.05/60)), size=0.2) + geom_rect(ymin=min(na.omit(dplot$signedminuslog_pval)), ymax=-log10(0.05/60), xmin = -Inf , xmax = Inf, alpha=0.002) + geom_bar(data=subset(dplot, minuslog_pval<(-log10(0.05/60))), aes(fill=log2_enrichment), fill="grey", stat="identity", colour='black')
    ## ggplot(dplot, aes(x = signedminuslog_pval, y = cell)) + geom_point(aes(fill=log2_enrichment), shape=21, colour='black', size=1.5) + facet_wrap(quantile~bin_number, nrow=1) + theme(plot.title=element_text(size=5), axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=6), panel.background = element_rect(fill="white"), legend.position="bottom", legend.justification="center", legend.title=element_text(size=5),legend.text=element_text(size=4), legend.key.size=unit(6,"mm") , panel.grid.major.y=element_line(size=0.3, colour='grey'), panel.grid.major.x=element_blank()) + labs(title=paste("eQTL enrichment in enhancer state clusters",i,sep=" "), y="cluster", x="Quantile bin") + scale_fill_gradientn(colours=c("lightblue","white","pink","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_vline(xintercept=c(min(na.omit(dplot$signedminuslog_pval)),-log10(0.05/60)), size=0.2) + geom_rect(xmin=min(na.omit(dplot$signedminuslog_pval)), xmax=-log10(0.05/60), ymin = -Inf , ymax = Inf, alpha=0.002)

})
dev.off()

## ## # Plot enrichment heatmap for Quintile bins
## denhancerCluster_quintile <- subset(denhancerCluster, quantile=="binQuintile")
## pdf(paste(args[2],"enhancerCluster_quintiles.pdf", sep="."), height=8, width=5)
## lapply(sort(unique(denhancerCluster_quintile$eqtl_set)), function(i){
##     dplot <- subset(denhancerCluster_quintile, eqtl_set==i)
##     ## dplotsignificant <- subset(dplot, minuslog <- pval>=2.79)
##     colVector <- c(min(dplot$log2_enrichment), 0, 0.7*(max(dplot$log2_enrichment)), max(dplot$log2_enrichment))
##     print(c(round(colVector,2)))
##     ggplot(dplot, aes(x = as.factor(bin_number), y = cell)) + geom_tile(data=dplot, aes(fill=log2_enrichment), colour='black') + facet_wrap(~quantile, nrow=1, scales="free_x") + theme(axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="black"), panel.grid=element_blank()) + labs(title=paste("eQTL enrichment in enhancer state clusters",i,sep=" "), y="cluster", x="Quantile bin") + scale_fill_gradientn(colours=c("lightblue","white","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_tile(data=subset(dplot, minuslog_pval<3.08), aes(fill=log2_enrichment),fill="grey") # geom_tile(subset(dplot, minuslog_pval<3.08), aes(x=as.factor(bin_number), y=cell), fill="grey")
## })
## dev.off()

## # Heatmap plots
## # Plot enrichment heatmap for fulleQTLset
## denhancerCluster_fullset <- subset(denhancerCluster, quantile=="fulleQTLset")
## pdf(paste(args[2],"enhancerCluster_fullset.pdf", sep="."), height=6, width=3)
## lapply(sort(unique(denhancerCluster_fullset$eqtl_set)), function(i){
##     dplot <- subset(denhancerCluster_fullset, eqtl_set==i)
##     ## dplotsignificant <- subset(dplot, minuslog <- pval>=2.79)
##     colVector <- c(min(dplot$log2_enrichment), 0, 0.7*(max(dplot$log2_enrichment)), max(dplot$log2_enrichment))
##     print(c(round(colVector,2)))
##     ggplot(dplot, aes(x = as.factor(bin_number), y = cell)) + geom_tile(data=dplot, aes(fill=log2_enrichment), colour='black') + facet_wrap(~quantile, nrow=1, scales="free_x") + theme(plot.title=element_text(size=5), axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=6), panel.background = element_rect(fill="black"), panel.grid=element_blank(), legend.position="bottom", legend.justification="center", legend.title=element_text(size=5),legend.text=element_text(size=4), legend.key.size=unit(4,"mm")) + labs(title=paste("eQTL enrichment in enhancer state clusters",i,sep=" "), y="cluster", x="Quantile bin") + scale_fill_gradientn(colours=c("lightblue","white","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_tile(data=subset(dplot, minuslog_pval<3.08), aes(fill=log2_enrichment),fill="grey")
## })
## dev.off()

## ## # Plot enrichment heatmap for Quintile bins
## denhancerCluster_quintile <- subset(denhancerCluster, quantile=="binQuintile")
## pdf(paste(args[2],"enhancerCluster_quintiles.pdf", sep="."), height=8, width=5)
## lapply(sort(unique(denhancerCluster_quintile$eqtl_set)), function(i){
##     dplot <- subset(denhancerCluster_quintile, eqtl_set==i)
##     ## dplotsignificant <- subset(dplot, minuslog <- pval>=2.79)
##     colVector <- c(min(dplot$log2_enrichment), 0, 0.7*(max(dplot$log2_enrichment)), max(dplot$log2_enrichment))
##     print(c(round(colVector,2)))
##     ggplot(dplot, aes(x = as.factor(bin_number), y = cell)) + geom_tile(data=dplot, aes(fill=log2_enrichment), colour='black') + facet_wrap(~quantile, nrow=1, scales="free_x") + theme(axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="black"), panel.grid=element_blank()) + labs(title=paste("eQTL enrichment in enhancer state clusters",i,sep=" "), y="cluster", x="Quantile bin") + scale_fill_gradientn(colours=c("lightblue","white","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_tile(data=subset(dplot, minuslog_pval<3.08), aes(fill=log2_enrichment),fill="grey") # geom_tile(subset(dplot, minuslog_pval<3.08), aes(x=as.factor(bin_number), y=cell), fill="grey")
## })
## dev.off()
