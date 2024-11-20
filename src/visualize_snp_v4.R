#!/usr/bin/env Rscript
#
# FT; last modification November 20th 2024
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
library(gridExtra) # multi graph on the same figure
library(cowplot)
library(stringr)
library(jsonlite) # to read json
b_verbose <- FALSE

args <- commandArgs(TRUE) # all arguments are character types
b_log_scale <- TRUE # default value, modified bar args[7]

density <- try( read.table(args[1], h=T, sep = "\t") ) # read the dataframe
if(inherits(density,"try-error"))
  density <- NULL

contig_limits <- try( read.table(args[2], h=F, sep = "\t") ) # read the dataframe
if(inherits(contig_limits,"try-error")){
  contig_limits <- NULL
}else
{
  contig_limits <- contig_limits$V1
}
  
contig_names <- try( read.table(args[3], h=F, sep = "\t") ) # read the dataframe
if(inherits(contig_limits,"try-error")){
  contig_names <- NULL
}else
{
  contig_names <- contig_names$V1

  # if contig name is \d+, add 'contig ' in front of it 
  replacement_contigname <- function(x){
    replaced = gsub("^(\\d+)$","contig \\1",x)            # replace :...., by ,
    return( replaced )
  }
  contig_names = lapply(contig_names, replacement_contigname)
  
  print(contig_names)
}  

threshold = as.numeric(args[4]) # define threshold

# added 2024 02 27
json_genes = read_json(args[5]) # json file with genes limits in "genes"->"name" [start,end]
if( b_verbose ){
  print("json genes:")
  str(json_genes$genes)
}

if( b_verbose ){
  print("json proteins:")
  str(json_genes$proteins)
}

outfile = args[6]

# added 2024 11 20: option to know if we want log10 scale for covdepth (1, true for the default; 0, false if not)
b_log_scale_int = as.numeric(args[7])
if( b_log_scale_int == 1)
{
  b_log_scale = TRUE
}else
{
  b_log_scale = FALSE
}


if( length(args) > 7 )
{
  b_covdepth <- TRUE
  # to prepare coverage depth graph above variant/annotation graph
  coverage_depth <- try( read.table(args[8], h=F, sep = "\t", blank.lines.skip=TRUE), silent=TRUE ) # read the dataframe
  if(inherits(coverage_depth,"try-error"))
    coverage_depth <- NULL
    # b_covdepth <- FALSE
    print("No coverage data found, disable coverage depth display in final graph")
}else
{ 
  b_covdepth <- FALSE 
}

width <- unit(21, "cm")

t1 = as.character(threshold) # define threshold as a character
# t = paste("Nucleotide Variation - threshold = ", t1, sep = ' ')

if(b_covdepth){

  if(! is.null(coverage_depth)){ # means covdepth provided, we prepare var for graph
    maxi = max(coverage_depth$V2)
    mini = min(coverage_depth$V2)
    med = median(coverage_depth$V2)
    m = paste("Coverage - median_coverage =", med, "[", mini, ":", maxi, "]", sep = " ")
    options(repr.plot.width = 5, repr.plot.height =1)
    if(b_log_scale)
    {
      cd = ggplot(coverage_depth, aes(x = V1, y = V2)) + geom_line(color = "black", linewidth = 0.5) + scale_x_continuous(trans='log10') + coord_trans(y="log10")
    }else
    {
      # when no log scale
      cd = ggplot(coverage_depth, aes(x = V1, y = V2)) + geom_line(color = "black", linewidth = 0.5)
    }
    #  + labs(list(title = m, x = "Base Position", y = "Number of Reads")
    cd1 = cd + labs(title = m) # add graph title
    cd2 = cd1 + xlab("") # add x axe title
    if(b_log_scale)
    {
      cd3 = cd2 + ylab("Number of reads (log)") # add axes and graph titles  
    }else
    {
      cd3 = cd2 + ylab("Number of reads") # add axes and graph titles
    } 
  }
  else{
    m = "Coverage - median_coverage = 0 [-inf:inf]"
    options(repr.plot.width = 5, repr.plot.height =1)
    cd = ggplot() + geom_line(color = "black", linewidth = 0.5)
    #  + labs(list(title = m, x = "Base Position", y = "Number of Reads")
    cd1 = cd + labs(title = m) # add graph title
    cd2 = cd1 + xlab("") # add x axe title
    cd3 = cd2 + ylab("Number of reads") # add axes and graph titles 
  }
}

options(repr.plot.width = 5, repr.plot.height =5) 
if(is.null(density)){ # means no snp found
    p = ggplot() # plot initialization
}else
{
    p = ggplot(density) # plot initialization
}

# defined a wide color palet manually to ensure better visibility, the other colors defined later will not be so distinctive
genecols=c('#99CC00','#CC9900','#FFCC33','#FF9900','#FF6600','#FF3300','#CC3300','#990033','#FF3366','#FF9999','#FFCCFF','#CC99CC','#996699','#993399','#660066','#990066','#660099','#CC00CC','#6600CC','#9966FF','#6633FF','#330099','#0033CC','#3366FF','#0033FF','#00CCFF','#006699','#00CC99','#009966','#33FFCC','#33CC99','#66FFCC','#33CCCC','#99FFFF','#66CC99','#339966','#99CC99','#99FF99','#66CC66','#66FF66','#669933','#336600','#00FF00','#66CC33','#CCFF66','#CCCC99','#FFFFCC','#FFCC99','#FFCCCC','#FF99CC','#CC99FF','#9999CC','#CCCCFF','#333333','#999999','#666666','#CCCCCC','#000000','#FF0000','#33D033','#CCD066','#CC00FF','#9900FF','#0099FF','#003333','#0099CC','#00FF99','#CCFF00')

# ,'#','#','#','#','#','#','#','#','#','#','#','#',
# '#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#',
# '#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#',
# '#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#','#')

# check color redundancy
genecols_unique = sort(unique(genecols))

if( b_verbose ){
  print("max_col_number_unique:")
  print(length(genecols_unique))
}

if( length(genecols) != length(genecols_unique) ){
    genecols = sort(genecols)
    for( icol in 1:length(genecols_unique) ){
      if( b_verbose ){
        print("genecols_unique de ")
        print(icol)
        print(":")
        print(genecols_unique[icol])
      }
      if( genecols_unique[icol] != genecols[icol] ){
        redund_vals = genecols[icol]
    	  print("Color found several times:")
    	  print(redund_vals)
    	  stop()
	}
    }
    print("first values identical, added values in genecols:")
    for( icol in length(genecols_unique)+1:length(genecols) ){
    	   print("Color found several times:")
    	   print(genecols[icol])
    }
    stop()
}

# add missing colors of rainbow 4000 to genecols, to be sure to always have enough colors (mimmivirus =~ 3000 genes)
complete_col_set=rainbow(n=10000)
indexes_to_rm = match(genecols, complete_col_set)
for(i in indexes_to_rm){
      complete_col_set[ -i ]
      }
genecols = append(genecols, complete_col_set)


# # check how many color we can manage
# print("max_col_number:")
# print(length(genecols))


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
gene_id_labels = unique(density$gene_id)
density$gene_id_num = paste(  sprintf("%02d", match(density$gene_id, gene_id_labels) ), ":", density$gene_id) 

# if( b_verbose ){
#   print("protein_id_labels BEFORE removing duplicates:")
#   print(density$protein_id)
# }

protein_id_labels = unique(density$protein_id)

# if( b_verbose ){
#   print("protein_id_labels BEFORE shortening:")
#   print(protein_id_labels)
# }

# truncate long legend
replacement <- function(x){
  replaced = str_replace_all(x, ":[[:alnum:] \\-]+,",",")            # replace :...., by ,
  replaced = str_replace_all(replaced, ":[[:alnum:] \\-]+:",":")     # replace :....: by :
  replaced = str_replace_all(replaced, ":[[:alnum:] \\-]+$","")      # remove last :....
  replaced = str_replace_all(replaced, "\\[[[:alnum:] \\-]+\\]","")  # remove [...]
  replaced = str_replace_all(replaced, " (?:putative|growth)[[:alnum:] \\-]+","")  # remove useless annotations
  replaced = str_replace_all(replaced, "similar to","~")  # remove useless annotations
  replaced = str_replace_all(replaced, "polyprotein","polyprot")  # remove useless annotations
  replaced = str_replace_all(replaced, "membrane","memb.")  # remove useless annotations
  replaced = str_replace_all(replaced, "envelope","env.")  # remove useless annotations
  replaced = str_replace_all(replaced, "protein","prot.")  # remove useless annotations
  replaced = str_replace_all(replaced, "coronavirus","coro.")  # remove useless annotations
  replaced = str_replace_all(replaced, " product,",",")  # remove useless annotations
  replaced = gsub(" hydrophobic domain [0-9]+", "", replaced)  # remove useless annotations
  replaced = str_replace_all(replaced, "polymerase","polym.")  # remove useless annotations
  replaced = str_replace_all(replaced, "dependent","dep.")  # remove useless annotations
  replaced = str_replace_all(replaced, "frameshift","fr.shift")  # remove useless annotations
  replaced = gsub(" metal binding [A-Za-z/]+","", replaced)  # remove useless annotations
  return( replaced )
}


# protein_id_labels = lapply(protein_id_labels, FUN = function(x) str_replace_all(x, ":[A-Za-z0-9 ]+,",""))
protein_id_labels = lapply(protein_id_labels, replacement)
if( b_verbose ){
  print("protein_id_labels SHORTENED:")
  print(protein_id_labels)
}
density$protein_id = lapply(density$protein_id, replacement)
# if( b_verbose ){
#   print("protein_id shortened:")
#   print(density$protein_id)
# }

density$protein_id_num = paste(  sprintf("%03d", match(density$protein_id, protein_id_labels) ), ":",  density$protein_id) 

# protein_id_ordered = fct_reorder(protein_id_labels, density$protein_id_num)

# # TODO stem_loop
# stem_loops_labels = unique(density$stem_loops)
# density$stem_loops_num = paste(  sprintf("%03d", match(density$stem_loops, stem_loops_labels) ), density$stem_loops) 

# print("gene_id_num:")
# print(density$gene_id_num)

if( b_verbose ){
  print("gene_id_labels:")
  print(gene_id_labels)
}

# add the consensus points, and remove the legend "linewidth" for proteins
p1 = p + geom_point(aes(x = position, y = -0.1, colour = density$gene_id_num, shape = density$protein_id_num), size = 1, show.legend = F, shape = 15) 
# add consensus gene boxes 
p1bis = p1

# -------------------------------------------------------------------------------
# for gene boxes display
gnames = names(json_genes$genes)
ybase  = -0.04
yshift = -0.00
xshift = -0.03
# minimal length of gene to consider to display name in gene boxes in the graph (bellow the absissia)
min_glength = 600

nrange = 1:length(gnames)

xmi=lapply(nrange, FUN=function(x){ json_genes$genes[[x]][[1]] } )
xmi=as.numeric(unlist(xmi))

xma=lapply(nrange, FUN=function(x){ json_genes$genes[[x]][[2]] } )
xma=as.numeric(unlist(xma))

# average x to write text into
xmm=lapply(nrange, FUN=function(x){ (xma[x] + xmi[x]) / 2 + xshift } )
xmm=as.numeric(unlist(xmm))

ymi=lapply(nrange, FUN=function(x){ ybase * (x %% 2) + yshift  } )
ymi=as.numeric(unlist(ymi)) # c(-0.04, -0.04,-0.04, -0.04,-0.04, -0.04,-0.04, -0.04,-0.04, -0.04,-0.04, -0.04,-0.04)

yma=lapply(nrange, FUN=function(x){ ybase * (1 + (x %% 2)) + yshift } ) 
yma=as.numeric(unlist(yma)) # c(-0.02, -0.02,-0.02, -0.02,-0.02, -0.02,-0.02, -0.02,-0.02, -0.02,-0.02, -0.02,-0.02)

# average y to write text into
ymm=lapply(nrange, FUN=function(x){ (yma[x] + ymi[x]) / 2 + yshift } )
ymm=as.numeric(unlist(ymm))

gnames_label=lapply(nrange, FUN=function(x){ 
                                    glabel = ""
                                    if( (xma[x] - xmi[x]) > min_glength){
                                      glabel=gnames[x]
                                    }
                                    glabel
                              } 
                      )
# colrect=rep("black", length(gnames))
# filrect=rep("white", length(gnames))
# -------------------------------------------------------------------------------

df=data.frame(xmi,xma,ymi,yma) # ,colrect,filrect)

if( b_verbose ){
  print(gnames)
  print(gnames_label)
  print(json_genes$genes[[1]][[1]])
  class(json_genes$genes[[1]][[1]])
  print(nrange)
  print(xmi)
  print(xma)
  print(xmm)
  print(ymi)
  print(yma)
  print(ymm)
  # print(colrect)
  # print(filrect)

  class(xmi)
  class(ymi)
  # class(colrect)
}

p1bis = p1 + geom_rect(data=df, mapping=aes(xmin=xmi, 
                                xmax=xma, 
                                ymin=ymi, 
                                ymax=yma), 
                                # color=colrect, 
                                # fill=filrect), 
                        color="black",
                        fill="white",
                        show.legend = FALSE,
                        inherit.aes = FALSE) + geom_text(
                                                  data = df,
                                                  aes(x = xmm, y = ymm, label = gnames_label),
                                                  size = 3, check_overlap = TRUE # vjust = 0, hjust = 0,
       )

# needed to have more than 6 symbols for gene_id legend
p1ter = p1bis + scale_shape_manual(values=c(1:25,1:25))

# use color more easily distinguishable
p1quad = p1ter + scale_color_manual(values=genecols[1:length(gene_id_labels)])

p2 = p1ter + geom_point(aes(x = position, y = as.numeric(variant_percent), color = density$gene_id_num, shape = density$protein_id_num)) # add the variant points
# p2 = p1ter + geom_point(aes(x = position, y = as.numeric(variant_percent), color = density$gene_id_num, shape = density$gene_id)) # add the variant points

p3 = p2 + geom_line(aes(x=position, y=threshold)) # add the threshold line

p3bis = p3			   
if(! is.null(contig_limits)){ # means more than 1 contig
  for(i in contig_limits){
        p3bis = p3bis + geom_vline(xintercept=i,linetype="dotted") # add the contig limits (vertical dotted line)
	}
  # average x to write text into
  namerangei = 1:length(contig_limits)

  if( b_verbose ){
    print("contig_limits:")
    print(contig_limits)
    print("contig_limits de 1:")
    print(contig_limits[1])
    print("contig_limits de length de contig_limits:")
    print(contig_limits[length(contig_limits)])
    print("namerangei:")
    print(namerangei)
  }

  xmidname=lapply(namerangei, FUN=function(x){ as.integer((as.integer(contig_limits[x-1]) + as.integer(contig_limits[x])) / 2) } )
  xmidname=as.numeric(unlist(xmidname))
  
  if( b_verbose ){
    print("xmidname:")
    print(xmidname)
  }
}
if(! is.null(contig_names)){ # means more than 1 contig
  if( b_verbose ){
    print(contig_names)
  }
  # for(ni in 1:(length(contig_names))){
  #   n = contig_names[ni]
  #   npos = xmidname[ni]
  #   print("n:")
  #   print(n)
  #   print("npos:")
  #   print(npos)
  #   p3bis = p3bis + geom_text(aes(x = npos, y = 1.1, label = n, angle = 0)) # add name to the graph
	# }
  cy = rep(1.1, length(contig_names)) 
  dfc = data.frame(xmidname,cy,contig_names)
  p3bis = p3bis + geom_text(data=dfc, 
                            mapping=aes(
                              x=xmidname, 
                              y=cy, 
                              label=contig_names),  
                            color="black",
                            size = 3, 
                            check_overlap = TRUE
                          )

}

# # TODO stem_loops
# if(! is.null(stem_loops)){ # means at least one stem_loop
#     for(i in stem_loops){
#         p4bis = p4bis + geom_point(aes(x=position, y=-0.15, colour=density$stem_loops_num) # add the stem loops
# 	}
# }


# p4 = p3bis + labs(title = t) # add graph title
p5 = p3bis + xlab("Base Position") # add x axe title
p6 = p5 + ylab("Variant Frequency") # add axes and graph titles

p6bis = p6 +guides( shape = guide_legend(order=2, 
                                        direction="vertical", 
                                        title="proteins (order: protein positions)", 
                                        nrow=16
                                        ), 
                    color = guide_legend(order=1, 
                                          direction="vertical", 
                                          title="genes (order: gene positions)"),
                                          nrow=3
                                          )

# modify the legend text sizes end position
p7 = p6bis + theme( plot.title = element_text(hjust=0.5),
                    legend.position="bottom",
                    legend.text = element_text(size=rel(0.6)), 
                    legend.title=element_text(size=rel(0.8)),
                    
                    ) # modify the legend position 

p9 = p7 + geom_text(aes(x = position, y = variant_percent + 0.03, label = indice, angle = 0)) # add indice to the graph
#p10 = p9 + geom_text(x = 0, y = threshold + 0.01, label = t) # add the threshold text
#p11 = p10 + geom_line(aes(x = position, y = 0.5), color = "red") 
p10 = p9 + geom_line(aes(x = position, y = 0.5), color = "red")
minus_maxlength_over6 = - max(contig_limits) / 50

# give scale and breaks of the y axis
p11 = p10 + scale_y_continuous(limits=c(-0.1,1.1), labels = scales::percent, breaks=c(0.0, 0.07, 0.2, 0.4, 0.6, 0.8, 1.0))

if( b_covdepth ){
  # gives proportions for covdepth graph (cd3) and variant graph (p11) in the grid plot
  g <- plot_grid(cd3, p11, align = "v", nrow = 2, rel_heights = c(1/6, 5/6))
}else
{
  # gives proportions for covdepth graph (cd3) and variant graph (p11) in the grid plot
  g <- plot(p11)
}

ggsave(outfile, device = "png", plot = g, width = width, units="cm", height=29.7, dpi = 600) # save the graph
# ~ end of script ~
