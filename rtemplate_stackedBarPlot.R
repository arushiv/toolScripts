# Stacked bar plot for % coverage of tiles in rat only, mouse only and both.
library(ggplot2)
library(reshape2)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=T)

#human_chromosome        total_tiles     mouse_only_mappings     percent_mouse_only_mappings    rat_only_mappings       percent_rat_only_mappings      mappings_mouseAndrat      percent_mappings_mouseAndrat     unmapped_tiles  percent_unmapped_tiles
df <- data.frame(d$human_chromosome, d$percent_mouse_only_mappings, d$percent_rat_only_mappings, d$percent_mappings_mouseAndrat, d$percent_unmapped_tiles)

colnames(df) <- c("chr","mouse_only","rat_only","mouse_rat_both","unmapped")
#df

dfm <- melt(df, id.vars=c("chr"))
colnames(dfm) <- c("chr", "mapping_type", "percent_maps_of_total_tiles")
#dfm

pdf(args[2], height=10, width=8)
ggplot(dfm, aes(x = chr, y = percent_maps_of_total_tiles, fill = mapping_type)) +  geom_bar(stat = "identity") + theme(text = element_text(size=12), axis.text.x = element_text(angle=45), panel.background = element_rect(fill = 'white')) 
dev.off()