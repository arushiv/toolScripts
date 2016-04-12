library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(stats)
library(qvalue)
library(ggrepel)
args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)
## d1 <- read.table(gzfile(args[3]), header=F)

# Put TFs into classes, show select TFs for enrichment
#eqtl    eqtl_set        quantile        bin_number      feature cell    subfeature      overlap expected_overlap        pval    tissue  motif_class

## selectMotifs <- d1[,1]
d$log2_enrichment <- log2(d$overlap/d$expected_overlap)
d$signedminuslog_pval <- -(log10(d$pval))*((d$log2_enrichment/abs(d$log2_enrichment)))
d$minuslog_pval <- -(log10(d$pval))
d <- subset(d, overlap>=as.numeric(args[3]))
qobj <- qvalue(d$pval)
## d$qsignificant <- qobj$significant
d$qval <- qobj$qvalues
d[,c("pval","qval")]
d$signedminuslog_qval <- -(log10(d$qval))*((d$log2_enrichment/abs(d$log2_enrichment)))
d$minuslog_qval <- -(log10(d$qval))

datacFootprint <- subset(d, feature%in%c("atacFootprints","allMotifs"))
## datacFootprint_isletSignificant <- subset(datacFootprint, subfeature%in%c(as.vector(d1[,1])))
## datacFootprint$subfeature <- factor(datacFootprint$subfeature, levels=rev(c(d1[,1]))
## datacFootprint_isletSignificant


## sort(unique(datacFootprint$cell))
# Plot enrichment heatmap for fulleQTLset
## datacFootprint_isletSignificant_fullset <- subset(datacFootprint_isletSignificant, quantile=="fulleQTLset")
bonferroni_threshold <- -log10(0.05/1751)
fdr_threshold <- -log10(0.05)

pdf(paste(args[2],"footprints_isletSignificant_classByName.qval.pdf", sep="."), height=12, width=10)
lapply(sort(unique(datacFootprint$cell)), function(i){
    d1 <- subset(datacFootprint, cell==i)
    lapply(sort(unique(d1$quantile)), function(j){
        d2 <- subset(d1, quantile==j)
        lapply(sort(unique(d2$eqtl_set)), function(k){
            dplot <- subset(d2, eqtl_set==k)
            ## dplotsignificant <- subset(dplot, minuslog <- pval>=2.79)
            ## colVector1 <- ifelse(min(d$signedminuslog_pval)<=bonferroni_threshold, c(min(dplot$signedminuslog_pval), -(bonferroni_threshold), 0, bonferroni_threshold, 0.35*(max(dplot$signedminuslog_pval)), 0.7*(max(dplot$signedminuslog_pval)), max(dplot$signedminuslog_pval))
            ## colVector <- c(min(dplot$signedminuslog_pval), 0, 0.35*(max(dplot$signedminuslog_pval)), 0.7*(max(dplot$signedminuslog_pval)), max(dplot$signedminuslog_pval))
            colVector <- c(min(dplot$signedminuslog_qval), 0, 0.35*(max(dplot$signedminuslog_qval)), 0.7*(max(dplot$signedminuslog_qval)), max(dplot$signedminuslog_qval))
            ## colVector <- c(min(dplot$log2_enrichment), 0, 0.35*(max(dplot$log2_enrichment)), 0.7*(max(dplot$log2_enrichment)), max(dplot$log2_enrichment))
            print(c(round(colVector,2)))
            ## log2(fold enrich on x axis)
            ## ggplot(dplot, aes(x = log2_enrichment, y = motif_class)) + geom_point(aes(fill=signedminuslog_pval), shape=21, colour='black', size=2) + facet_wrap(bin_number~feature, nrow=1) + theme(axis.title.y=element_text(size=5), axis.title.x=element_text(size=5), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="white"), panel.grid.major.y=element_line(size=0.2, colour='black'), panel.grid.major.x=element_blank(), strip.text.x=element_text(size=5)) + labs(title=paste("eQTL enrichment in TF classes",i,j,k,sep=" "), y="Motif class", x="-log10(p value)") + scale_fill_gradientn(colours=c("lightblue","white","pink","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_point(data=subset(dplot, minuslog_pval<(bonferroni_threshold)), aes(fill=signedminuslog_pval), shape=21, fill="grey") + geom_text_repel(data=subset(dplot, minuslog_pval>=bonferroni_threshold), aes(label=subfeature), size=1, colour='black', segment.size=0.2, force=2, box.padding=unit(.3,'lines'), point.padding=unit(.15,'lines'))#+ geom_vline(xintercept=c(min(na.omit(dplot$signedminuslog_pval)),-log10(0.05/1751)), size=0.2) + geom_rect(xmin=min(na.omit(dplot$signedminuslog_pval)), xmax=-log10(0.05/1751), ymin = -Inf , ymax = Inf, alpha=0.02)# geom_tile(subset(dplot, minuslog_pval<3.08), aes(x=as.factor(bin_number), y=cell), fill="grey")

            ## log10(pval) on x axis, bonferroni correction
            ## ggplot(dplot, aes(x = signedminuslog_pval, y = motif_class)) + geom_point(aes(fill=log2_enrichment), shape=21, colour='black', size=2) + facet_wrap(bin_number~feature, nrow=1) + theme(axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="white"), panel.grid.major.y=element_line(size=0.2, colour='black'), panel.grid.major.x=element_blank()) + labs(title=paste("eQTL enrichment in TF classes",i,j,k,sep=" "), y="Motif class", x="-log10(p value)") + scale_fill_gradientn(colours=c("lightblue","white","pink","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_vline(xintercept=c(min(na.omit(dplot$signedminuslog_pval)),bonferroni_threshold), size=0.2) + geom_rect(xmin=min(na.omit(dplot$signedminuslog_pval)), xmax=bonferroni_threshold, ymin = -Inf , ymax = Inf, alpha=0.02) + geom_text_repel(data=subset(dplot, minuslog_pval>=bonferroni_threshold), aes(label=subfeature), size=1, colour='black', segment.size=0.2, force=2, box.padding=unit(.3,'lines'), point.padding=unit(.15,'lines')) # geom_tile(subset(dplot, minuslog_pval<3.08), aes(x=as.factor(bin_number), y=cell), fill="grey")

            ## log10(qvalue) on x axis after FDR 
            ## ggplot(dplot, aes(x = signedminuslog_qval, y = motif_class)) + geom_point(aes(fill=log2_enrichment), shape=21, colour='black', size=2) + facet_wrap(bin_number~feature, nrow=1) + theme(axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="white"), panel.grid.major.y=element_line(size=0.2, colour='black'), panel.grid.major.x=element_blank()) + labs(title=paste("eQTL enrichment in TF classes",i,j,k,sep=" "), y="Motif class", x="-log10(q value)") + scale_fill_gradientn(colours=c("lightblue","white","pink","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_vline(xintercept=c(min(na.omit(dplot$signedminuslog_qval)), fdr_threshold), size=0.2) + geom_rect(xmin=min(na.omit(dplot$signedminuslog_qval)), xmax=fdr_threshold, ymin = -Inf , ymax = Inf, alpha=0.02) + geom_text_repel(data=subset(dplot, minuslog_qval>=fdr_threshold), aes(label=subfeature), size=1, colour='black', segment.size=0.2, force=2, box.padding=unit(.3,'lines'), point.padding=unit(.15,'lines')) # geom_tile(subset(dplot, minuslog_pval<3.08), aes(x=as.factor(bin_number), y=cell), fill="grey")


            ## log2(fold enrich on x axis) after FDR
            ggplot(dplot, aes(x = log2_enrichment, y = motif_class)) + geom_point(aes(fill=signedminuslog_qval), shape=21, colour='black', size=2) + facet_wrap(bin_number~feature, nrow=1) + theme(axis.title.y=element_text(size=5), axis.title.x=element_text(size=5), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="white"), panel.grid.major.y=element_line(size=0.2, colour='black'), panel.grid.major.x=element_blank(), strip.text.x=element_text(size=5)) + labs(title=paste("eQTL enrichment in TF classes",i,j,k,sep=" "), y="Motif class", x="-log10(q value)") + scale_fill_gradientn(colours=c("lightblue","white","pink","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) + geom_point(data=subset(dplot, minuslog_qval<(fdr_threshold)), aes(fill=signedminuslog_qval), shape=21, fill="grey") + geom_text_repel(data=subset(dplot, minuslog_qval>=fdr_threshold), aes(label=subfeature), size=1, colour='black', segment.size=0.2, force=2, box.padding=unit(.3,'lines'), point.padding=unit(.15,'lines'))#+ geom_vline(xintercept=c(min(na.omit(dplot$signedminuslog_pval)),-log10(0.05/1751)), size=0.2) + geom_rect(xmin=min(na.omit(dplot$signedminuslog_pval)), xmax=-log10(0.05/1751), ymin = -Inf , ymax = Inf, alpha=0.02)# geom_tile(subset(dplot, minuslog_pval<3.08), aes(x=as.factor(bin_number), y=cell), fill="grey")

        })
    })
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
##     ggplot(dplot, aes(x = as.factor(bin_number), y = cell)) + geom_tile(data=dplot, aes(fill=log2_enrichment), colour='black') + facet_wrap(~quantile, nrow=1, scales="free_x") + theme(axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="black"), panel.grid=element_blank()) + labs(title=paste("eQTL enrichment in enhancer state clusters",i,sep=" "), y="cluster", x="Quantile bin") + scale_fill_gradientn(colours=c("lightblue","white","red","darkred"), values=rescale(c(colVector)), breaks=c(colVector), labels=c(round(colVector,2))) # geom_tile(subset(dplot, minuslog_pval<3.08), aes(x=as.factor(bin_number), y=cell), fill="grey")
## })
## dev.off()
