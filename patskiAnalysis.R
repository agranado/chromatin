#Analysis of Patski cell line ATAC data for epimemoir
#Friday, 12th October,
#Aim: identify long contiguous regions of open/closed chromatin for designing gRNA targets



library(rtracklayer)

library(IdeoViz)
library(GenomicRanges)
library(readxl)
#setwd("/home/alejandrog/MEGA/Caltech/epimemoir/chromatin") #GIT repository
patski= import.bed("patskiWT.allelicATAC.dScore.minCov_5.bedgraph.gz")


  #Distribution of region size
  x11()
  hist(log10(end(ranges(patski))-start(ranges(patski))),
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
patski.gr.unique = all.segmented.gr[[1]] #replace the object with the partitioned genome.
#let's take only the ATAC peaks from the X chromosome
patski.X.segmented = patski.gr.unique[ seqnames(patski.gr.unique)=="chrX" ]
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



#data frame can be filtered by value, in this case a high score reprensents the regions we want to keep:
score.threshold = 0.5
new.regions<-x[x$value>score.threshold,]
#Make a GenomicRanges object
patski.open.gr<-makeGRangesFromDataFrame(new.regions,ignore.strand = T,seqnames.field = "chrom",keep.extra.columns = T)


# # # # # # # # #

# # # # # # # # #
 # # # # # # # # #
### FIND GENES BASED ON REGIONS

# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf

source("mus.musculusX.R")
#order by score
coord.promoters.x = promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,upstream=2000,downstream=0)
coord.genes.x = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

#make a data frame from the promoters GR object (keeping all fields)
promoters.df<-data.frame(chrom=seqnames(coord.promoters.x),start=start(coord.promoters.x),
                          end=end(coord.promoters.x),strand =strand(coord.promoters.x) ,tx_id=coord.promoters.x$tx_id,
                          tx_name=coord.promoters.x$tx_name)

#alternative approach :
  #for each promoter, look for overlapping ATAC peaks /scores
  raw.patski = all.bed.files[[1]] #extract the patksi ATAC data unlisted (basically the original data with no duplicates)
  raw.patski.x = raw.patski[seqnames(raw.patski)=="chrX",] #only genes from chr X
  high.patski.x=raw.patski.x[raw.patski.x$score>0.45]
  summary(width(high.patski.x)) # what's the average length of the ATAC peaks in patski data

  atac.promoters = overlap_ATAC_ranges( promoters.df, raw.patski.x )
#this promoters have UC ID which is weird, we can make a conversion table to get gene symbols and mrnaID
#> mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D mm10 -e 'select kgID,mRNA,geneSymbol from kgXref' | sed "s/'/\'/;s/\t/,/g;s/^//;s/$//;s/\n//g" > mm10_geneSymbols.txt
#then read (needs a couple of manuall edits cause some entries have multiple names)
gene.ids=read.table("mm10_geneSymbols.txt",sep=",",header=T)
#load the genome:
library(BSgenome.Mmusculus.UCSC.mm10)
chrx<-unmasked(Mmusculus$chrX)
#get sequences for top promoters
promoter.seqs = Views(chrx,start = atac.promoters$start,end = atac.promoters$end)
library(stringi)
#add the sequences for all the promoters
atac.promoters$sequence = as.character(promoter.seqs)
#get reverse complement
for(i in 1:dim(atac.promoters)[1]){
  if(atac.promoters[i,]$strand=="-"){
    atac.promoters[i,]$sequence = rc(atac.promoters[i,]$sequence)
  }
}

# alternative
# >> tx_seqs1<-extractTranscriptSeqs(Mmusculus,TxDb.Mmusculus.UCSC.mm10.knownGene,use.names = T)
