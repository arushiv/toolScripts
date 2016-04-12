library(ggplot2)
library(reshape2)
library(scales)
args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)


                                        # colour assignments
                                        # Stretch enhancer labelling
                                        # If it lies in an ATAC-seq peak

## cell    chromatin_state trait snp_chr snp_start       snp_end SNP     Proxy   Distance        RSquared        DPrime  Arrays  Major   Minor   MAF     NObserved atac_peak
## d <- subset(d, trait=="T2D" & SNP=="rs1535500")

d <- subset(d, trait==args[3])
d$snpInfo <- paste(d$Proxy, d$snp_chr, d$snp_start, d$snp_end, sep="_")
d1 <- subset(d, chromatin_state!="stretchEnhancer" & chromatin_state!="typicalEnhancer" & chromatin_state!="intermediateEnhancer" & chromatin_state!="stretch4_2kbEnhancer" & chromatin_state!="stretch6kbEnhancer")
da <- subset(d1, atac_peak=="1")
d1$chromatin_state <- factor(d1$chromatin_state, levels=c("1_Active_TSS","2_Weak_TSS","3_Flanking_TSS","5_Strong_transcription","6_Weak_transcription","8_Genic_enhancer","9_Active_enhancer_1","10_Active_enhancer_2","11_Weak_enhancer","14_Bivalent_poised_TSS","16_Repressed_polycomb","17_Weak_repressed_polycomb","18_Quiescent_low_signal"))

d2 <- subset(d,chromatin_state=="stretchEnhancer" | chromatin_state=="typicalEnhancer" | chromatin_state=="intermediateEnhancer")

pdf(args[2], height=7, width=8)
lapply(sort(unique(d1$SNP)), function(i){
    ggplot(subset(d1, SNP==i), aes(x = snpInfo, y = cell)) + geom_tile(aes(fill=chromatin_state), colour='black') + facet_wrap(~SNP, scales="free_x", ncol=1) + theme(axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), text = element_text(size=5), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="black"), panel.grid=element_blank()) + scale_fill_manual(values=c("5_Strong_transcription"=rgb(0,128,0,maxColorValue=255),"6_Weak_transcription"=rgb(0,100,0,maxColorValue=255),"8_Genic_enhancer"=rgb(194,225,5,maxColorValue=255),"11_Weak_enhancer"=rgb(255,255,0,maxColorValue=255),"9_Active_enhancer_1"=rgb(255,195,77,maxColorValue=255),"10_Active_enhancer_2"=rgb(255,195,77,maxColorValue=255),"1_Active_TSS"=rgb(255,0,0,maxColorValue=255),"3_Flanking_TSS"=rgb(255,69,0,maxColorValue=255),"2_Weak_TSS"=rgb(255,69,0,maxColorValue=255),"14_Bivalent_poised_TSS"=rgb(205,92,92,maxColorValue=255),"18_Quiescent_low_signal"=rgb(255,255,255,maxColorValue=255),"17_Weak_repressed_polycomb"=rgb(192,192,192,maxColorValue=255),"16_Repressed_polycomb"=rgb(128,128,128,maxColorValue=255))) + geom_point(data=subset(d2, SNP==i), aes(shape=chromatin_state), fill='white', colour='black', size=1) + geom_point(data=subset(da, SNP==i), shape=8, size=2, colour='black')
})
dev.off()

## Original single plot ##
## n <- length(unique(d1$SNP))
## pdf(args[2], height=(7*n), width=8)
##  ## for (i in 1:length(vector)){           
## ggplot(d1, aes(x = snpInfo, y = cell)) + geom_tile(aes(fill=chromatin_state), colour='black') + facet_wrap(~SNP, scales="free_x", ncol=1) + theme(axis.title.y=element_text(size=7), axis.title.x=element_text(size=7), text = element_text(size=5), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y=element_text(size=7), panel.background = element_rect(fill="black"), panel.grid=element_blank()) + scale_fill_manual(values=c("5_Strong_transcription"=rgb(0,128,0,maxColorValue=255),"6_Weak_transcription"=rgb(0,100,0,maxColorValue=255),"8_Genic_enhancer"=rgb(194,225,5,maxColorValue=255),"11_Weak_enhancer"=rgb(255,255,0,maxColorValue=255),"9_Active_enhancer_1"=rgb(255,195,77,maxColorValue=255),"10_Active_enhancer_2"=rgb(255,195,77,maxColorValue=255),"1_Active_TSS"=rgb(255,0,0,maxColorValue=255),"3_Flanking_TSS"=rgb(255,69,0,maxColorValue=255),"2_Weak_TSS"=rgb(255,69,0,maxColorValue=255),"14_Bivalent_poised_TSS"=rgb(205,92,92,maxColorValue=255),"18_Quiescent_low_signal"=rgb(255,255,255,maxColorValue=255),"17_Weak_repressed_polycomb"=rgb(192,192,192,maxColorValue=255),"16_Repressed_polycomb"=rgb(128,128,128,maxColorValue=255))) + geom_point(data=d2, aes(shape=chromatin_state), size=1) + geom_point(data=da, shape=8, size=2, colour='black', fill='white')
## dev.off()


## pdf(args[3], height=6, width=6)
## 	for (i in 1:length(vector)){
## 	maxE <- max(subset(d1,footprint %in% c(vector[i]))$enrich, na.rm=T)
## 	minE <- min(subset(d1,footprint %in% c(vector[i]))$enrich, na.rm=T)
## 	color_scale <- round(c(minE, 0, 0.33*maxE, .66*maxE, maxE), digits=2)
## 	print(ggplot(subset(d1, footprint %in% c(vector[i])), aes(x = chromatin_state, y = state_tissue)) + facet_grid(type~region) + geom_tile(aes(fill=enrich)) + scale_fill_gradientn(colours=c("lightblue","white","pink", "red", "darkred"), values=rescale(color_scale), breaks=color_scale, labels=color_scale, name="enrichment=log2(real_stat/mean_null)") + theme(text = element_text(size=8), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=6), axis.text.y = element_text(size=4), panel.background = element_rect(fill="white")) + labs(title=paste("Enrichments for chromatin states overlaps with ATAC-seq",vector[i], sep=" "), x="Chromatin States", y="Cell/Tissue type") + scale_x_discrete(limits=x_order))
## 	}
## 	dev.off()
