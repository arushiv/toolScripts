library(ggplot2)
library(reshape2)
library(scales)
library(ggrepel)
args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=T)
d1 <- read.table(gzfile(args[2]), header=T)

# ensemblId       islet_fpkm      entropy isletSpecificity        entropySpec     biotype gene_name       iesi    decile
# gene    adipose adrenal brain   breast  colon   heart   islet   kidney  liver   lung    lymph_node      ovary   prostate        skeletal_muscle testes  thyroid white_blood_cells       entropy isletSpecificity        entropySpec     biotype class   name

diesi <- subset(d, select = c("ensemblId","gene_name","decile","iesi"))
drpkm <- subset(d1, select = c("gene", "name", "adipose", "adrenal", "brain", "breast", "colon", "heart", "islet", "kidney", "liver", "lung", "lymph_node", "ovary", "prostate", "skeletal_muscle", "testes", "thyroid", "white_blood_cells"))
## names(drpkm)
colnames(diesi) <- c("gene","name","islet_decile","iesi")
## names(diesi)

df <- merge(diesi, drpkm, by=c("gene","name"))
## names(df)
d2 <- melt(df, id=c("gene","name","islet_decile","iesi"))
colnames(d2) <- c("gene","name","islet_decile","iesi","tissue","rpkm")

# Get quintile from decile
l1 <- c(1:10)
l2 <- c(1,1,2,2,3,3,4,4,5,5)
l <- cbind(l1,l2)

d2$islet_quintile <- unlist(lapply(d2$islet_decile, function(j){
    apply(l, 1, function(i){
        if(j==i[1]){
            return(i[2])
        }
    })
}), use.names=FALSE)

## d2[,c("islet_decile","islet_quintile")]
tomatch <- c("RFX")

d2 <- subset(d2, tissue %in% c("islet", "adipose", "brain", "heart", "liver", "skeletal_muscle"))
d2$tissue <- factor(d2$tissue, levels=c("adipose", "brain", "heart", "liver", "skeletal_muscle", "islet"))
dplot <- subset(d2, rpkm!=0)
## dplot[, c("islet_decile","islet_quintile")] 
## Plot scatter for iESI vs RPKM
pdf(args[3], height=4, width=10)
ggplot(dplot, aes(x=rpkm, y=iesi)) + geom_point(aes(fill=as.factor(islet_quintile), colour=as.factor(islet_quintile)), shape=21, size=1, alpha=0.6) + facet_wrap(~tissue, nrow=1, scales="free_x") + theme(legend.position="bottom", text = element_text(size=7), axis.text.x = element_text(vjust=1, hjust=0.5, size=4), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="iESI vs RPKM in tissues across iESI quintile bins", x="Gene RPKM in respective tissue", y="islet expression specificity index (iESI)") + scale_x_log10() + geom_label_repel(data=subset(dplot, rpkm>=1000), aes(label=name), size=1, segment.size=0.2, force=2, box.padding=unit(.1,'lines'), point.padding=unit(.15,'lines'), fill='white', color='black') + geom_label_repel(data=subset(subset(dplot, tissue=="islet")[grep(tomatch, subset(dplot, tissue=="islet")$name),], islet_decile=="10"), aes(label=name), size=1, segment.size=0.2, force=2, box.padding=unit(.1,'lines'), point.padding=unit(.15,'lines'), fill='white', color='red') + scale_fill_manual(name="Quintile bin", values=c("blue","lightblue","yellow","orange","red")) + scale_color_manual(name="Quintile bin", values=c("blue","lightblue","yellow","orange","red"))
dev.off()


## Plot iESI denisity

## dplot <- subset(dplot, tissue=="islet")
## lmin <- unlist(lapply(sort(unique(dplot$islet_quintile)), function(i){
##     return(min(subset(dplot, islet_quintile==i)$iesi))
## }), use.names=FALSE)
## lmax <- unlist(lapply(sort(unique(dplot$islet_quintile)), function(i){
##     return(max(subset(dplot, islet_quintile==i)$iesi))
## }), use.names=FALSE)
## lcol <- c("blue","lightblue","yellow","orange","red")
## lmat <- cbind(lmin,lmax, lcol)
## dens <- density(dplot$iesi)
## dd <- with(dens,data.frame(x,y))


## ## p <- ggplot(subset(dplot, tissue=="islet"), aes(iesi)) + geom_density() + theme(panel.background = element_rect(fill = 'white', colour = 'black')) + coord_flip()
## p <- qplot(x,y,data=dd,geom="line") + theme(panel.background = element_rect(fill = 'white', colour = 'black')) + coord_flip()
## p <- p + unlist(apply(lmat, 1, function(i){
##     ## with(dens, polygon(x=c(i[1], i[1]:i[2], i[2]), y=c(0, y[i[1]:i[2]], 0), col=i[3]))
##     geom_ribbon(data=subset(dd, x>=i[1] & x<=i[2]), aes(ymax=y, ymin=0) , fill=i[3], colour=NA, alpha=0.5)
## }), use.names=FALSE)
## pdf(args[3], height=4, width=2)
## plot(p)
## dev.off()


## d1 <- d1[, !(names(d1) %in% c("genename","gene"))]
## colnames(d1) <- c("binQuintile","tissue","rpkm")
## d2 <- aggregate(rpkm ~ binQuintile + tissue, data=d1, mean, na.rm=TRUE)
## colnames(d2) <- c("binQuintile","tissue","mean_rpkm")

## dliver <- subset(d1, tissue=="liver", binQuintile=="1")
## dislet <- subset(d1, tissue=="islets", binQuintile=="1")
## dtest <- subset(d1, tissue=="liver" | tissue=="islets" | tissue=="testes" | tissue=="adrenal" | tissue=="skeletal_muscle")
## dtest1 <- subset(d1, tissue=="liver" | tissue=="islets" | tissue=="testes")
## max(dliver$rpkm)
## max(dislet$rpkm)
## minE <- min(d2$mean_rpkm)
## maxE <- max(d2$mean_rpkm)

## pdf(args[2], height=13, width=11)
## ## ggplot(d2, aes(x=tissue, y=binQuintile)) + geom_tile(aes(fill=mean_rpkm)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="Mean RPKM across iESI quintile bins across tissues") + scale_fill_gradientn(colours=c("blue","lightblue","yellow","red"), limits=c(minE,maxE))
## ## Boxplot:
## ## ggplot(d1, aes(factor(tissue), rpkm)) + geom_boxplot(aes(fill=as.factor(binQuintile))) + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="RPKM for genes in iESI quintile bins across tissues") + scale_y_log10()
## ggplot(d1, aes(y=binQuintile, x=rpkm)) + geom_point(aes(fill=tissue), alpha=0.5, shape=21) + facet_wrap(~tissue, ncol=1) + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="RPKM for genes in iESI quintile bins across tissues") + scale_x_log10()
## #+ scale_fill_gradientn(colours=c("blue","lightblue","yellow","red"), limits=c(minE,maxE))
## # ggplot(dboxplot, aes(factor(region_fraction), ks_D))
## ## ggplot(dtest1, aes(rpkm)) + geom_density(aes(fill=as.factor(binQuintile)), alpha=0.5) + facet_wrap(~tissue) + theme(strip.text.y=element_text(size=5), axis.text.y=element_text(size=5), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + scale_x_log10()
## ## ggplot(dtest, aes(rpkm)) + geom_density(aes(fill=tissue), alpha=0.5) + facet_grid(~binQuintile) + theme(strip.text.y=element_text(size=5), axis.text.y=element_text(size=5), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + scale_x_log10()
## ## ggplot(d1, aes(rpkm)) + geom_density(aes(fill=tissue), alpha=0.5) + facet_grid(~binQuintile) + theme(strip.text.y=element_text(size=5), axis.text.y=element_text(size=5), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + scale_x_log10()

## dev.off()


