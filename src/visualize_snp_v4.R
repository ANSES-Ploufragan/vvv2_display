#!/usr/bin/env Rscript
#
# AF; last modification February 16th 2019
#
# Description:  This script has been written in order to generate a graph
#               representing the viral genome and the possible variants inside
#               the viral population.
#               This script has been written in order to be executed as the last
#               step of a pipeline named Me.
#               It uses ggplot2 in order to draw this graph.
#
################################################################################
# ~ start of script ~

library(ggplot2) # import the ggplot2 library
args <- commandArgs(TRUE) # all arguments are character types

density <- try( read.table(args[1], h=T, sep = "\t") ) # read the dataframe
if(inherits(density,"try-error"))
  density <- NULL

threshold = as.numeric(args[2]) # define threshold

t1 = as.character(threshold) # define threshold as a character
t = paste("Nucleotide Variation - threshold = ", t1, sep = ' ')

if(is.null(density)){ # means no snp found
  p = ggplot() # plot initialization
}else{
  p = ggplot(density) # plot initialization
}

p1 = p + geom_point(aes(x = position, y = 0, color = gene_id), size = 1, show.legend = F, shape = 15) # add the consensus points, and remove the legend "size"  
p2 = p1 + geom_point(aes(x = position, y = as.numeric(variant_percent), color = gene_id)) # add the variant points
p3 = p2 + geom_line(aes(x=position, y=threshold)) # add the threshold line
p4 = p3 + labs(title = t) # add graph title
p5 = p4 + xlab("Base Position") # add x axe title
p6 = p5 + ylab("Variant Frequency") # add axes and graph titles
p7 = p6 + theme(legend.position = "bottom") # modify the legend position
p8 = p7 + ylim(-0.06,1.2) # modify the scale
p9 = p8 + geom_text(aes(x = position, y = variant_percent + 0.01, label = indice, angle = 0)) # add indice to the grapĥ
#p10 = p9 + geom_text(x = 0, y = threshold + 0.01, label = t) # add the threshold text
#p11 = p10 + geom_line(aes(x = position, y = 0.5), color = "red") 
p10 = p9 + geom_line(aes(x = position, y = 0.5), color = "red")
p11 = p10 + theme(plot.title = element_text(hjust=0.5))


ggsave(args[3], device = "png", plot = last_plot(), width = 20, dpi = 600) # save the graph

# ~ end of script ~
