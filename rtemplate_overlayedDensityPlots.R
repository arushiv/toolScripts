library(ggplot2)
library(reshape2)
library(scales)
library(ggrepel)
library(gridExtra)
library(gapminder)

args <- commandArgs(TRUE)
dcorr <- read.table(gzfile(args[1]), header=T)
dpval <- read.table(gzfile(args[2]), header=T)

                                        # tissue tf1 tf2 ...
                                        # tissue tf1 tf2 ...

bonferroni_threshold <- -log10(0.05/60)
tomatch <- c("RFX", "MA0509", "MA0510", "MA0600")

pdf(args[3], height=7, width=6)
lapply(sort(unique(dcorr$tissue)), function(i){
    d1 <- cbind(as.character(as.vector(names(dcorr))), as.numeric(as.vector(subset(dcorr, tissue==i))))
    d2 <- cbind(as.character(as.vector(names(dpval))),as.numeric(as.vector(subset(dpval, tissue==i))))
    colnames(d1) <- c("tf","correlation")
    colnames(d2) <- c("tf","pval")
    d <- merge(d1,d2,by="tf")
    d <- subset(d, tf!="tissue")
    print(typeof(d$pval))
    d$minuslog_pval <- -log10(as.numeric(as.character(d$pval)))
    subset(d, minuslog_pval>=5.9)

    p1 <- ggplot(d, aes(x=as.numeric(as.character(correlation)), y=minuslog_pval)) + geom_point(stat="identity", fill="grey", shape=21, size=1, alpha=0.6) + xlim(-1,1) + ylim(0,6.005) +
          geom_point(data=subset(d, minuslog_pval>=bonferroni_threshold), stat="identity", fill="red", shape=21, size=1, alpha=0.3) +
          theme(plot.title=element_text(size=5), axis.text.x = element_text(size=6), panel.background = element_rect(fill = 'white', colour = 'black'), axis.title.x=element_blank()) +
          labs(title=paste("tissue=",i,sep=" "), y="-log10(p value)") +
          geom_label_repel(data=subset(d, minuslog_pval>=5.9), aes(label=tf), size=1, segment.size=0.2, force=2, box.padding=unit(.1,'lines'), point.padding=unit(.15,'lines'), fill='white', color='red') +
          geom_label_repel(data=d[grep(paste(tomatch, collapse="|"), d$tf),], aes(label=tf), size=1, segment.size=0.2, force=2, box.padding=unit(.1,'lines'), point.padding=unit(.15,'lines'), fill='white', color='black') +
          geom_vline(aes(xintercept=0), colour='black', size=0.1) +
          geom_hline(aes(yintercept=bonferroni_threshold), colour='red', size=0.2)

    p2 <- ggplot(d, aes(x=as.numeric(as.character(correlation)))) + geom_density() + xlim(-1,1) +
          theme(panel.background=element_rect(fill='white', colour='black'),  axis.text.x = element_text(size=6)) +
          labs(x="Correlation")

    p3 <- ggplot(d, aes(x=minuslog_pval)) + geom_density() + coord_flip() + xlim(0,6) +
          theme(panel.background=element_rect(fill='white', colour='black'), axis.title=element_blank(), plot.title=element_text(size=5), axis.text.x = element_text(size=6)) +
          labs(title=paste("tissue=",i,sep=" ")) +
          geom_vline(aes(xintercept=bonferroni_threshold), colour='red', size=0.2)

    p4 <- ggplot(d, aes(x=as.numeric(as.character(correlation)))) + geom_density(colour='white') + xlim(-1,1) +
        theme(panel.background=element_rect(fill='white', colour='white'), axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank())
    
    grid.arrange(p1,p3,p2,p4, nrow=2, heights=c(.80,.20), widths=c(.80,.20))    
})

dev.off()    


## Original plot for islets:

## d1 <- cbind(as.character(as.vector(names(dcorr))), as.numeric(as.vector(subset(dcorr, tissue=="Islets"))))
## d2 <- cbind(as.character(as.vector(names(dpval))),as.numeric(as.vector(subset(dpval, tissue=="Islets"))))
## colnames(d1) <- c("tf","correlation")
## colnames(d2) <- c("tf","pval")

## d <- merge(d1,d2,by="tf")
## d <- subset(d, tf!="tissue")
## typeof(d$pval)
## d$minuslog_pval <- -log10(as.numeric(as.character(d$pval)))
## subset(d, minuslog_pval>=5.9)
## bonferroni_threshold <- -log10(0.05/60)
## max(d$minuslog_pval)
## tomatch <- c("RFX", "MA0509", "MA0510", "MA0600")

## pdf(args[3], height=7, width=6)
## ggplot(d, aes(x=as.numeric(as.character(correlation)), y=minuslog_pval)) + geom_point(stat="identity", colour=NULL, fill="grey", shape=21, size=1, alpha=0.6) + xlim(-1,1) + ylim(0,6.2) +
##     geom_point(data=subset(d, minuslog_pval>=bonferroni_threshold), stat="identity", colour=NULL, fill="red", shape=21, size=1, alpha=0.6) +
##     theme(axis.text.x = element_text(size=6), panel.background = element_rect(fill = 'white', colour = 'black')) +
##     labs(x="Correlation", y="-log10(p value)") +
##     geom_label_repel(data=subset(d, minuslog_pval>=5.9), aes(label=tf), size=1, segment.size=0.2, force=2, box.padding=unit(.1,'lines'), point.padding=unit(.15,'lines'), fill='white', color='red') +
##     geom_label_repel(data=d[grep(paste(tomatch, collapse="|"), d$tf),], aes(label=tf), size=1, segment.size=0.2, force=2, box.padding=unit(.1,'lines'), point.padding=unit(.15,'lines'), fill='white', color='black') +
##     geom_vline(aes(xintercept=0), colour='black', size=0.1) +
##     geom_hline(aes(yintercept=bonferroni_threshold), colour='red', size=0.2)
## dev.off()


