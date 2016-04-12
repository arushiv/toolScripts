library(ggplot2)
library(reshape2)
library(scales)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

#Cell_1  Region_1        Cell_2  Region_2        region_fraction	subRegion_fraction        replicate       mean_null       mean    normalized_bootstrap_sd real_stat       scaled_real_stat        z_score pval    ks_D    ks_pval

mean_D <- aggregate( ks_D ~ Cell_1 + Region_1 + Cell_2 + Region_2 + region_fraction + subRegion_fraction, data = d, mean)
sd_D <- aggregate( ks_D ~ Cell_1 + Region_1 + Cell_2 + Region_2 + region_fraction + subRegion_fraction, data = d, sd)
d1 <- merge(mean_D, sd_D, by=c("Cell_1","Region_1","Cell_2","Region_2","region_fraction", "subRegion_fraction"))
vector <- levels(unique(d1$Region_2))	

pdf(args[2], height=10, width=10)
for (i in 1:length(vector)){
dsubset <- subset(d1, Region_2 %in% c(vector[i]))
minE <- min(dsubset$ks_D.x, na.rm=T)	
maxE <- max(dsubset$ks_D.x, na.rm=T)
color_scale <- c(maxE, minE+.5*(maxE-minE) , minE+0.33*(maxE-minE), minE)
print(color_scale)
dmin <- subset(dsubset, ks_D.x == minE)
print(ggplot(subset(d1, Region_2 %in% c(vector[i])), aes(x = as.character(region_fraction), y = as.character(subRegion_fraction))) + geom_tile(aes(fill=ks_D.x)) + geom_point(data=subset(dsubset, ks_D.x == minE), color = "green") + scale_fill_gradientn(colours=c("white", "yellow", "pink", "darkred"), values=rescale(color_scale), name="Average min D statistic for 15 GSC runs") + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=20), axis.text.y = element_text(size=20), panel.background = element_rect(fill="white")) +  scale_colour_hue(l=10, c=250) + labs(title=paste("ks test D after GSC for GM12878 stretch enhancer overlap with",vector[i], sep=" "), x="Region_fraction", y="Subregion_fraction"))
	}
dev.off()
	
pdf(args[3], height=6, width=8)
dboxplot <- subset(d, Region_2 == "broadDomains" & subRegion_fraction == "0.1")
ggplot(dboxplot, aes(factor(region_fraction), ks_D)) + geom_boxplot() + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="k-s D vs r for GM12878 stretch enhancer overlap with broadDomains -s=0.1, n=15")  
dev.off()







	
# ggplot(d, aes(x=regional_overlap, y=ks_D)) + geom_line() + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="k-s statistic D vs r for GSC runs (x in log2 scale)") + scale_x_continuous(trans=log2_trans())
# ggplot(d, aes(x=r, y=D)) + geom_line() + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=1, hjust=1), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="k-s statistic D vs r for GSC runs (x in log10 scale)") + scale_x_log10(limits=c(0.0001,0.7))
# ggplot(d, aes(x=r, y=D)) + geom_line() + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=1, hjust=1,  size=5), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="k-s statistic D vs r for GSC runs") + scale_x_continuous(breaks=seq(0, 1, by=0.01))
# ggplot(d, aes(x=r, y=pval)) + geom_line() + theme(text = element_text(size=10), axis.text.x = element_text(angle=45, vjust=0.5, hjust=0.5), panel.background = element_rect(fill = 'white', colour = 'black')) + labs(title="k-s statistic pval vs r for GSC runs") + scale_x_continuous(breaks=seq(0, 1, by = .025))
# dev.off()

