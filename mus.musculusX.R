
#from https://www.biostars.org/p/167818/
 # if if is not installed:  : >biocLite("Mus.musculus")
 library(Mus.musculus)


 #set the X chromosome as the only active for analysis
 seqlevels(TxDb.Mmusculus.UCSC.mm10.knownGene)<-"chrX"

# this object is automatically created by the library
 coord.genes.x = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
 coord.promoters.x = promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,upstream=2000,downstream=400)

 gene.length.x = width(ranges(coord.genes.x))
