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
colnames(diesi) <- c("gene","name","islet_decile","iesi")

df <- merge(diesi, drpkm, by=c("gene","name"))
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

## Plot iESI denisity
dplot <- subset(dplot, tissue=="islet")
lmin <- unlist(lapply(sort(unique(dplot$islet_quintile)), function(i){
    return(min(subset(dplot, islet_quintile==i)$iesi))
}), use.names=FALSE)
lmax <- unlist(lapply(sort(unique(dplot$islet_quintile)), function(i){
    return(max(subset(dplot, islet_quintile==i)$iesi))
}), use.names=FALSE)
lcol <- c("blue","lightblue","yellow","orange","red")
lmat <- cbind(lmin,lmax, lcol)
dens <- density(dplot$iesi)
dd <- with(dens,data.frame(x,y))

p <- qplot(x,y,data=dd,geom="line") + theme(panel.background = element_rect(fill = 'white', colour = 'black')) + coord_flip()
p <- p + unlist(apply(lmat, 1, function(i){
    geom_ribbon(data=subset(dd, x>=i[1] & x<=i[2]), aes(ymax=y, ymin=0) , fill=i[3], colour=NA, alpha=0.5)
}), use.names=FALSE)

pdf(args[3], height=4, width=2)
plot(p)
dev.off() 
