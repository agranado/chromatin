#Analysis of Patski cell line ATAC data for epimemoir
#Friday, 12th October,
#Aim: identify long contiguous regions of open/closed chromatin for designing gRNA targets



library(rtracklayer)

library(IdeoViz)
library(GenomicRanges)
library(readxl)
setwd("/home/alejandrog/MEGA/Caltech/epimemoir/chromatin") #GIT repository
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

#There are some repetead entries in the data
patski.gr.unique = unique(patski.gr)
#Let's partition the genome into equal length long regions
#We can use a custom function
#input argument is a list of bed objects
all.bed.files  = list()
all.bed.files[[1]] = patski.gr.unique #only element in the list

#Load the function :
source("grFunctions.R")

#return a genomic ranges object
#We still need to modify this

#partition the genome in windows to look for density of atac peaks
window.size=5000
all.segmented.gr<-segment_genome(all.bed.files,window.size)

#let's take only the ATAC peaks from the X chromosome
patski.X = patski.gr.unique[ seqnames(patski.gr.unique)=="chrX" ]
#WE can se how many sites we for a give score
# There are 1605 peaks in the X chromose of which 204 have score 0.5 (maximum)
#Sushi library
library(Sushi)
x<-data.frame(chrom = seqnames(patski.X.segmented),start = start(patski.X.segmented),
            end = end(patski.X.segmented),value = patski.X.segmented$score)

chrom = "chrX";chromstart = min(x$start);chromend=max(x$end)
plotBedgraph(x,chrom = "chrX",chromstart = min(x$start),chromend = max(x$end))
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("ATAC d-score",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
