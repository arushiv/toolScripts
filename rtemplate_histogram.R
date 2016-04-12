library(ggplot2)


args <- commandArgs(TRUE)
df <- read.table(args[1], header=F)
colnames(df) <- c("mouse_lengths", "rat_lengths")

png('mouse_lengths.png')
hist(df$mouse_lengths, xlab="lengths of hg19 tiles mappings to mouse", main="Distribution of lengths of mm9 sequences that map to 200bp hg19 tiles")
dev.off()

vector <- levels(unique(df$cell))
pdf(args[3], height=12, width=12)
for (i in 1:length(vector)){
maxE <- max(subset(df,cell %in% c(vector[i]))$lengths, na.rm=T)
minE <- min(subset(df,cell %in% c(vector[i]))$lengths, na.rm=T)
print(ggplot(subset(df,cell %in% c(vector[i])), aes(lengths, fill=stretchEnh_model)) + geom_density(alpha=0.5, binwidth=200) + labs(title=paste("Stretch Enhancer length distribution for Parker et. al. and 13 State chromatin state models for", vector[i], sep=" "))+ theme(strip.text.x = element_text(size = 11), axis.text.x = element_text(angle=90, vjust=0.5, hjust=1) + panel.background = element_rect(fill = 'white')) + scale_x_continuous(limits = c(minE, maxE)))
}
dev.off()
