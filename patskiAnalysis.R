#Analysis of Patski cell line ATAC data for epimemoir
#Friday, 12th October, 
#Aim: identify long contiguous regions of open/closed chromatin for designing gRNA targets 



library(rtracklayer)

library(IdeoViz)
library(GenomicRanges)
library(readxl)

patski= import.bed("patskiWT.allelicATAC.dScore.minCov_5.bedgraph.gz")


#Distribution of region size 
x11()
hist(log(end(ranges(patski))-start(ranges(patski))),
     main="Region size, all chromosomes",xlab="Region size (log bp)")

#Most regions seem to be less than 1kb
#median = 669

#reformat as data.frame to include numeric score values 
patski.df = data.frame(chr = seqnames(patski), start = start(ranges(patski)), 
                       end= end(ranges(patski)), score = as.numeric( patski$name) )
#remake genomic ranges object 
patski.gr=makeGRangesFromDataFrame(patski.df,ignore.strand = T,
                                    seqnames.field = "chr",keep.extra.columns = T)

#Let's partition the genome into equal length long regions 
#We can use a custom function
#input argument is a list of bed objects
all.bed.files  = list()
all.bed.files[[1]] = patski.gr #only element in the list

#Load the function : 
source("../mouseATAC/grFunctions.R")

#return a genomic ranges object
window.size=5000
all.segmented.gr<-segment_genome(all.bed.files,window.size)


