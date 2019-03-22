
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


#same function as before but when the GR object has actual numerical scores
grToDataFrameScore<-function(chronis.atac){
  closed.val = 0
  open.val=1
  chronis.score = rep(1,length(chronis.atac))
  #make data frame with same info (helps to process the data)
  celln1<-data.frame( seqnames=seqnames(chronis.atac),start=start(chronis.atac)-1,
                      end=end(chronis.atac),
                      score=chronis.atac$score,type=c(rep("Intergenic",length(chronis.score))))
  #Rename the colums such that my old script works
  names(celln1)<-c("chr","start","end","score","type")
  #define the levels 0->closed; 1->open
  #levels(celln1$score)<-c(0,1)
  return(celln1)
}


#modified for X patski cell line analysis:
#this is part of an independent GIT repository for anlaysis of patski cell line

#major changes:
#   Add new function grToDataFrameScore that actually includes the numerical score provided in the original data
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
      celln1.list[[i]]<-grToDataFrameScore(all.atac[[i]])

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
          all.scores= array()

          for(c in 1:length(all.chr)){
                chr=all.chr[c]
                #how many windows of 10kb do we have in this chromosome
                segmented.chr.size=round((min.ends[c] - max.starts[c])/window.size)
                #the new array
                segmented.chromosome =array(0,segmented.chr.size)
                new.starts = array(0,segmented.chr.size)
                new.ends= array(0,segmented.chr.size)
                new.scores = array(0,segmented.chr.size)

                #get all the start / end regions from this chromosome
                st.chr=celln1$start[celln1$chr==chr]
                en.chr = celln1$end[celln1$chr==chr]
                score.chr = celln1$score[celln1$chr==chr]

                init.region =max.starts[c]
                end.region = max.starts[c] + window.size +1

                for(s in 1:length(segmented.chromosome)){
                    #Here en.chr includes all the endpoints of the regions (peaks) in the original data.
                    # segmented.chromosome will COUNT the total number of peaks within each NEW region of the
                    # segmented data ( by window.size)

                    #if the data had score we could include them here:
                    segmented.chromosome[s] = sum(en.chr<= end.region & st.chr>= init.region)
                    #new values: # # # # # #
                    new.starts[s] = init.region
                    new.ends[s] = end.region
                    #since we are mergin peaks, lets's
                    #take the scores for all peaks in this regions and mutiplied by their values
                    #then add them as weighted sum
                    if(segmented.chromosome[s]>0){
                      new.scores[s] = sum(score.chr[en.chr<= end.region & st.chr>= init.region])
                    }
                    # end new values # # # # # #

                    #update indexes:
                    init.region = init.region + window.size +1
                    end.region = init.region + window.size
                }
                seq.names=rep(chr,length(segmented.chromosome))

                all.seq.names = c(all.seq.names,seq.names)
                all.new.starts = c(all.new.starts,new.starts)
                all.new.ends = c(all.new.ends, new.ends)
                all.segmented = c(all.segmented,segmented.chromosome)
                all.scores = c(all.scores, new.scores)
          }
          segmented.dataframe=data.frame(chr=all.seq.names,start=all.new.starts,end=all.new.ends,score=all.scores,npeaks = all.segmented)
          segmented.dataframe=segmented.dataframe[!is.na(segmented.dataframe$start),]
          segmented.gr=makeGRangesFromDataFrame(segmented.dataframe,ignore.strand=T,seqnames.field="chr",keep.extra.columns=T)

          all.segmented.gr[[i]]=segmented.gr
     }



    return(all.segmented.gr)



}
#reverse complement for ORF on the negative strand
rc <- function(nucSeq){
  #return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))
  return(stri_reverse(nucSeq))
}

#this function takes annotated promoters from UCSD browser and then looks for peaks (ATAC in this case) within the promoter regions.
#target.regions is a dataframe with start end attributes
#ATAC peaks is a genomic region with start, end and score
overlap_ATAC_ranges<-function(target.regions,ATACpeaks){
  raw.patski.x = ATACpeaks
  promoters.df = target.regions
  promoter.score = array(0,dim(promoters.df)[1])
  promoter.peaks = array(0,dim(promoters.df)[1])
  promoter.mean = array(0,dim(promoters.df)[1])
  promoter.peaks.location = array(0,dim(promoters.df)[1])

  for(i in 1:dim(promoters.df)[1]){
    promoter.peaks[i] = sum( start(raw.patski.x)>=promoters.df[i,]$start & end(raw.patski.x)<=promoters.df[i,]$end)
    promoter.score[i] = sum(raw.patski.x$score[start(raw.patski.x)>=promoters.df[i,]$start & end(raw.patski.x)<=promoters.df[i,]$end])
    promoter.mean[i] =  mean(raw.patski.x$score[start(raw.patski.x)>=promoters.df[i,]$start & end(raw.patski.x)<=promoters.df[i,]$end])


    #FIND location of ATAC peaks and convert the location relative to the start of the sequnece
    #ALSO reverse the sequence if it is on the negative strand (and fix the location relative to the END!)
    s1=start(raw.patski.x)[start(raw.patski.x)>=promoters.df[i,]$start & end(raw.patski.x)<=promoters.df[i,]$end]
    s2=end(raw.patski.x)[start(raw.patski.x)>=promoters.df[i,]$start & end(raw.patski.x)<=promoters.df[i,]$end]
    location =""
    if(length(s1)>0){
        ### FIX location of the peaks
        if(promoters.df[i,]$strand=="+"){
        s1= s1-promoters.df[i,]$start
        s2 = s2-promoters.df[i,]$start
      }else { #promoter is in the reverse strand
        s1_ = promoters.df[i,]$end - s2
        s2 = promoters.df[i,]$end - s1
        s1=  s1_
      }


        j=1;while(j <= length(s1)){
          location=paste(location,toString(s1[j]),",",toString(s2[j]),";",sep="")
          j=j+1
        }
    }



    promoter.peaks.location[i] = location

  }
  promoters.df$sum.score = promoter.score
  promoters.df$n.peaks = promoter.peaks
  promoters.df$mean.score  = promoter.mean
  promoters.df$peaks.location = promoter.peaks.location

#NEW data.frame for all promoters with 2 new fields:
#  score:  (sum of individual scores for peaks within the promoter)
#  peaks:  number of peaks ATAC peaks found within the region
  atac.promoters<-promoters.df[promoters.df$peaks>0,]
  atac.promoters<-atac.promoters[base::order(atac.promoters$sum.score,decreasing=T),]
  return(atac.promoters)
}


cluster_chromosome<-function(all.segmented.gr,chr){

  chr.matrix=c()
  for(l in 1:length(all.segmented.gr)){
        this.chr=all.segmented.gr[[l]][ seqnames(all.segmented.gr[[l]])==chr  ]$score
        chr.matrix = cbind(chr.matrix,this.chr)
    }

}
