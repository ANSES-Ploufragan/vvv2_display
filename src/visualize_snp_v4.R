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

# install.packages("viridis")
# library("viridis")
library(ggplot2) # import the ggplot2 library
args <- commandArgs(TRUE) # all arguments are character types

density <- try( read.table(args[1], h=T, sep = "\t") ) # read the dataframe
if(inherits(density,"try-error"))
  density <- NULL

contig_limits <- try( read.table(args[2], h=F, sep = "\t") ) # read the dataframe
if(inherits(contig_limits,"try-error"))
  contig_limits <- NULL
  
threshold = as.numeric(args[3]) # define threshold

t1 = as.character(threshold) # define threshold as a character
t = paste("Nucleotide Variation - threshold = ", t1, sep = ' ')

if(is.null(density)){ # means no snp found
  p = ggplot() # plot initialization
}else{
  p = ggplot(density) # plot initialization
}

# defined a wide color palet manually to ensure better visibility, the other colors defined later will not be so distinctive
protcols=c('#99CC00','#CC9900','#FFCC33','#FF9900','#FF6600','#FF3300','#CC3300','#990033','#FF3366','#FF9999','#FFCCFF','#CC99CC','#996699','#993399','#660066','#990066','#660099','#CC00CC','#6600CC','#9966FF','#6633FF','#330099','#0033CC','#3366FF','#0033FF','#00CCFF','#006699','#00CC99','#009966','#33FFCC','#33CC99','#66FFCC','#33CCCC','#99FFFF','#66CC99','#339966','#99CC99','#99FF99','#66CC66','#66FF66','#669933','#336600','#00FF00','#66CC33','#CCFF66','#CCCC99','#FFFFCC','#FFCC99','#FFCCCC','#FF99CC','#CC99FF','#9999CC','#CCCCFF','#333333','#999999','#666666','#CCCCCC','#000000','#FF0000','#33D033','#CCD066','#CC00FF','#9900FF','#0099FF','#003333','#0099CC','#00FF99','#CCFF00')

# ,'#','#','#','#','#','#','#','#','#','#','#','#',
# '#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#',
# '#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#',
# '#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#')

# check color redundancy
protcols_unique = sort(unique(protcols))
print("max_col_number_unique:")
print(length(protcols_unique))

if( length(protcols) != length(protcols_unique) ){
    protcols = sort(protcols)
    for( icol in 1:length(protcols_unique) ){
       print("protcols_unique de ")
       print(icol)
       print(":")
       print(protcols_unique[icol])
       if( protcols_unique[icol] != protcols[icol] ){
           redund_vals = protcols[icol]
    	   print("Color found several times:")
    	   print(redund_vals)
    	   stop()
	}
    }
    print("first values identical, added values in protcols:")
    for( icol in length(protcols_unique)+1:length(protcols) ){
    	   print("Color found several times:")
    	   print(protcols[icol])
    }
    stop()
}

# add missing colors of rainbow 4000 to protcols, to be sure to always have enough colors (mimmivirus =~ 3000 genes)
complete_col_set=rainbow(n=10000)
indexes_to_rm = match(protcols, complete_col_set)
for(i in indexes_to_rm){
      complete_col_set[ -i ]
      }
protcols = append(protcols, complete_col_set)


# check how many color we can manage
print("max_col_number:")
print(length(protcols))


# print("genes:")
# print(density$gene_id)
# print("\n\n")

# print("proteins:")
# print(density$protein_id)
# print("\n\n")

## --------------------------------
## TODO stem_loops RNAs
# print("rnas:")
# print(density$rnas)
# print("\n\n")

# print("stem_loops:")
# print(density$stem_loops)
# print("\n\n")
## --------------------------------

# if(! is.null(contig_limits)){ # means 1 contig only
#      print("contig_limits:")
#      print(contig_limits)
#      print("\n\n")
# }

# needed to have color labels NOT ORDERED and to choose colors
protein_id_labels = unique(density$protein_id)
density$protein_id_num = paste(  sprintf("%03d", match(density$protein_id, protein_id_labels) ), density$protein_id) 

# # TODO stem_loop
# stem_loops_labels = unique(density$stem_loops)
# density$stem_loops_num = paste(  sprintf("%03d", match(density$stem_loops, stem_loops_labels) ), density$stem_loops) 

# print("protein_id_num:")
# print(density$protein_id_num)

# print("protein_id_labels:")
# print(protein_id_labels)

p1 = p + geom_point(aes(x = position, y = -0.05, colour = density$protein_id_num, shape = density$gene_id), size = 1, show.legend = F, shape = 15) # add the consensus points, and remove the legend "size" for proteins

# needed to have more than 6 symbols for gene_id legend
p1bis = p1 + scale_shape_manual(values=c(1:25))

# use color more easily distinguishable
p1ter = p1bis + scale_color_manual(values=protcols[1:length(protein_id_labels)])

p2 = p1ter + geom_point(aes(x = position, y = as.numeric(variant_percent), color = density$protein_id_num, shape = density$gene_id)) # add the variant points

p3 = p2 + geom_line(aes(x=position, y=threshold)) # add the threshold line

p3bis = p3			   
if(! is.null(contig_limits)){ # means more than 1 contig
    for(i in contig_limits){
        p3bis = p3bis + geom_vline(xintercept=i,linetype="dotted") # add the contig limits (vertical dotted line)
	}
}

# # TODO stem_loops
# if(! is.null(stem_loops)){ # means at least one stem_loop
#     for(i in stem_loops){
#         p4bis = p4bis + geom_point(aes(x=position, y=-0.15, colour=density$stem_loops_num) # add the stem loops
# 	}
# }


p4 = p3bis + labs(title = t) # add graph title
p5 = p4 + xlab("Base Position") # add x axe title
p6 = p5 + ylab("Variant Frequency") # add axes and graph titles

p6bis = p6 +guides(shape = guide_legend(order=1, direction="horizontal", title="gene"),
                   color = guide_legend(order=2, direction="vertical", title="protein"))
		   
# p7 = p6 + theme(legend.position = "bottom") # modify the legend position
p8 = p6bis + ylim(-0.06,1.2) # modify the scale
p9 = p8 + geom_text(aes(x = position, y = variant_percent + 0.03, label = indice, angle = 0)) # add indice to the grapĥ
#p10 = p9 + geom_text(x = 0, y = threshold + 0.01, label = t) # add the threshold text
#p11 = p10 + geom_line(aes(x = position, y = 0.5), color = "red") 
p10 = p9 + geom_line(aes(x = position, y = 0.5), color = "red")
p11 = p10 + theme(plot.title = element_text(hjust=0.5))


ggsave(args[4], device = "png", plot = last_plot(), width = 20, dpi = 600) # save the graph

# ~ end of script ~
