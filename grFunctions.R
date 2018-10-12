
grToDataFrame<-function(chronis.atac){
  closed.val = 0
  open.val=1
  chronis.score = rep(1,length(chronis.atac))
  #make data frame with same info (helps to process the data)
  celln1<-data.frame( seqnames=seqnames(chronis.atac),start=start(chronis.atac)-1,
                      end=end(chronis.atac),
                      scores=c(rep(open.val,length(chronis.atac))),type=c(rep("Intergenic",length(chronis.score))))
  #Rename the colums such that my old script works 
  names(celln1)<-c("chr","start","end","score","type")
  #define the levels 0->closed; 1->open
  levels(celln1$score)<-c(0,1)
  return(celln1)
}


segment_genome<-function(all.atac,window.size){
  #all.atac is a list of GR objects
 
  celln1.list=list()
  all.chr = levels(seqnames(all.atac[[1]])) 
   

    starts.chr = matrix(0,length(all.atac),length(all.chr))
    end.chr= matrix(0,length(all.atac),length(all.chr))
    
    #Here we will generate a matrix of data.frames, each one from an ATAC seq dataset
    #all the samples in the list will be clustered in the future, so we need to normalize their length
    for(i in 1:length(all.atac)){
      #convert each ATAC dataset to a dataframe
      celln1.list[[i]]<-grToDataFrame(all.atac[[i]])
      
      celln1.gr=all.atac[[i]]
      all.chr=levels(seqnames(celln1.gr))
      #lets find out the range of each chromosome for each data set
      for(c in 1:length(all.chr)){
        celln1 = celln1.list[[i]]
        chr = all.chr[c]
        starts.chr[i,c] = min(celln1$start[celln1$chr==chr  ]   )
        end.chr[i,c] = max(celln1$end[celln1$chr == chr  ]  )
      }
    }
    
    #look for the max value acrros the starting points: 
    #we need to start all samples at the region in the chromosome, 
    #so we need to start from the one with the highest value
    max.starts = apply(starts.chr,2,max)
    min.ends = apply(end.chr,2,min)
    
    
    #now that we have the range of the samples we can partition the genome in contigous 10kb windows
    
    
    all.segmented.gr = list()
    #lets go through all the datasets in the list,
     for(i in 1:length(all.atac)){   
          all.chr=levels(seqnames(all.atac[[i]]))
          celln1 = celln1.list[[i]]
          all.seq.names = array()
          all.new.starts= array()
          all.new.ends = array()
          all.segmented = array()
          
          for(c in 1:length(all.chr)){
                chr=all.chr[c]
                #how many windows of 10kb do we have in this chromosome
                segmented.chr.size=round((min.ends[c] - max.starts[c])/window.size)
                #the new array 
                segmented.chromosome =array(0,segmented.chr.size)
                new.starts = array(0,segmented.chr.size)
                new.ends= array(0,segmented.chr.size)
     
                #get all the start / end regions from this chromosome
                st.chr=celln1$start[celln1$chr==chr]
                en.chr = celln1$end[celln1$chr==chr]
                
                init.region =max.starts[c]
                end.region = max.starts[c] + window.size +1
                
                for(s in 1:length(segmented.chromosome)){
    
                    segmented.chromosome[s] = sum(en.chr<= end.region & st.chr>= init.region)
                    new.starts[s] = init.region
                    new.ends[s] = end.region
                      
                    init.region = init.region + window.size +1
                    end.region = init.region + window.size
                }
                seq.names=rep(chr,length(segmented.chromosome))
                
                all.seq.names = c(all.seq.names,seq.names)
                all.new.starts = c(all.new.starts,new.starts)
                all.new.ends = c(all.new.ends, new.ends)
                all.segmented = c(all.segmented,segmented.chromosome)
 
          }
          segmented.dataframe=data.frame(chr=all.seq.names,start=all.new.starts,end=all.new.ends,score=all.segmented)
          segmented.dataframe=segmented.dataframe[!is.na(segmented.dataframe$start),]
          segmented.gr=makeGRangesFromDataFrame(segmented.dataframe,ignore.strand=T,seqnames.field="chr",keep.extra.columns=T)
          
          all.segmented.gr[[i]]=segmented.gr
     }  
    
    
    
    return(all.segmented.gr)
    
    
    
}


cluster_chromosome<-function(all.segmented.gr,chr){
  
  chr.matrix=c()
  for(l in 1:length(all.segmented.gr)){
        this.chr=all.segmented.gr[[l]][ seqnames(all.segmented.gr[[l]])==chr  ]$score
        chr.matrix = cbind(chr.matrix,this.chr)
    }
  
}