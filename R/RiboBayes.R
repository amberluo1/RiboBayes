library(tidyverse)
library(ribor)
library(psych)
library(wavethresh)
library(rstatix)
library(sgr)
library(parallel)
library(Metrics)
library(reshape2)
library(edgeR)

pkm.ribo <- Ribo("/Users/amberluo/Downloads/pkm.ribo", rename = rename_default )
mm.ribo <- Ribo("/Users/amberluo/Downloads/all (1).ribo", rename = rename_default )
all.ribo <- Ribo("/Users/amberluo/Downloads/all.ribo", rename = rename_default )
# CC -> Added range lower/upper for consistency
lengths=lengths_aggregated(pkm.ribo)
mmlengths=lengths_aggregated(mm.ribo)
shifts=get_pshifts(ribo, graph=FALSE)

rc_CDS <- get_region_counts(ribo.object    = ribo,
                            range.lower = lengths[[1]],
                            range.upper = lengths[[length(lengths)]],
                            tidy       = TRUE,
                            alias      = TRUE,
                            transcript = FALSE,
                            normalize=TRUE,
                            region     = "CDS",
                            compact    = FALSE)
region_lengths <- get_internal_region_lengths(ribo.object = ribo, alias = TRUE)
cds=region_lengths%>%select(transcript, CDS)
rc_CDS=rc_CDS%>%left_join(cds)%>%mutate(count=count/CDS)

rc_CDS_w = dcast(rc_CDS[,-5], transcript ~ experiment)
# A value of 500 gives 283 transcripts; 1500 gives 106. 
high_exp = rowSums( cpm(rc_CDS_w[,-1]) > 369) > 1
sum(high_exp)

#highexp=rc_CDS%>%filter(count>0.4125, experiment=="20201209-shRNA-PKM-KD-1-ribo")
transcripts=rc_CDS_w$transcript[high_exp]
experiments=get_experiments(ribo)

rc_CDS <- get_region_counts(ribo.object    = mm.ribo,
                            range.lower = mmlengths[1],
                            range.upper = mmlengths[length(mmlengths)],
                            tidy       = TRUE,
                            alias      = TRUE,
                            transcript = FALSE,
                            normalize=TRUE,
                            region     = "CDS",
                            compact    = FALSE)
#retrieves the 350 most expressed transcripts
# CC removed deprecated function; updated selection of highexpression genes to include all experiments
region_lengths <- get_internal_region_lengths(ribo.object = mm.ribo, alias = TRUE)
cds=region_lengths%>%dplyr::select(transcript, CDS)
rc_CDS=rc_CDS%>%left_join(cds)%>%mutate(count=count/CDS)

rc_CDS_w = dcast(rc_CDS[,-5], transcript ~ experiment)
# A value of 500 gives 283 transcripts; 1500 gives 106. 
high_exp = rowSums( cpm(rc_CDS_w[,-c(1, 8:11)]) > 416 ) > 2
sum(high_exp)
mmtranscripts=rc_CDS_w$transcript[high_exp]
#cleaned
wavetransform=function(cov){
  experiments=unique(cov$experiment)
  max=max(cov$position)
  power=floor(log(max, 2))
  noise=c(rep(0.001,nrow(cov)))
  noise=jitter(noise, factor = 1, amount=NULL)
  wavecov=cov%>%arrange(experiment)%>%ungroup()%>%mutate(noises=noise)
  transformed=c()
  for(i in 1:length(experiments)){
    subcov=wavecov%>%filter(experiment %in% experiments[i])%>%mutate(counts=count+noises)
    noisecounts=subcov$counts
    noisecounts1=noisecounts[1:(2^power)]
    noisecounts2=noisecounts[(max-(2^power-1)):max]
    blocks.thr1 <- wavethresh::BAYES.THR(noisecounts1, plotfn=FALSE, filter.number=1, 
                                         family = "DaubExPhase", alpha=0, beta=0)
    blocks.thr2 <- wavethresh::BAYES.THR(noisecounts2, plotfn=FALSE, filter.number=1, 
                                         family = "DaubExPhase", alpha=0, beta=0)
    if(all(blocks.thr1==blocks.thr2)){
      test=blocks.thr1
    }else{
    test=c(blocks.thr1, blocks.thr2[(2^(power+1)+1-max):(2^power)])
    }
    transformed=c(transformed, test)
  }
  wavecov[,4]=transformed
  return(wavecov)
}
#clean
makeset=function(zcov, zcov_filtered){
  if(nrow(zcov_filtered)==0){
    return(0)
  }
  nrows=nrow(zcov_filtered)
  max=max(zcov_filtered$position)
  positions=unique(zcov_filtered$position)
  count=length(positions)
  mypositions=zcov%>%filter(position %in% positions)
  set=dcast(mypositions, position ~ experiment,value.var="score")
  numexperiments=ncol(set)-1
  set=set%>%mutate(avg1=rowMeans(set[,c(2:(numexperiments/2+1))]), avg2=rowMeans(set[,c((numexperiments/2+2):(numexperiments+1))]))
  m=1
  for(i in 1:nrow(set)){
    if(set[[m,(2+numexperiments)]]<7 & set[[m,(3+numexperiments)]]<7){
      set=set[-m,]
      m=m-1
    }
    m=m+1
  }
  if(nrow(set)==0){
    return(0)
  }
  positions=set$position
  mypositions=zcov%>%group_by(experiment)%>%mutate(originalz=(count-mean(count))/sd(count))%>%filter(position %in% positions)
  set=dcast(mypositions, position ~ experiment, value.var="originalz")
  return(set)
}
#same as makeset, but replaces the scores with 1 for pause site and 0 for not pause site
makekappaset=function(zcov, zcov_filtered, threshold){
  nrows=nrow(zcov_filtered)
  
  zcov_filtered=zcov_filtered%>%arrange(desc(position))
  max=zcov_filtered[[1,2]]
  
  vector3 <- c(1:max)
  count=0
  for(i in 1:max){
    if(vector3[i]==0){
      count=count+1
    }
  }
  vector4=c(1:count)
  for(i in 1:nrows){
    a=zcov_filtered[[i, 2]]
    vector3[a]=0
  }
  count=1
  for(i in 1:max){
    if(vector3[i]==0){
      vector4[count]=i
      count=count+1
    }
  }
  count=count-1
  set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
                 "control1"=1:count, "control2"=1:count, "control3"=1:count)
  #vector4 records all nucleotide sites in the GAPDH transcript with score > 30 in at
  #least one of the experiments
  #constructing the skeleton for set
  for(i in 1:count){
    set[i, 1]=vector4[i]
  }
  
  for(i in 1:count){
    for(j in 2:7){
      set[i,j]=0
    }
  }
  
  zcov3=zcov%>%filter(position%in%vector4)
  count=1
  zcov3=zcov3%>%arrange(position)
  
  rows=nrow(set)
  #fills in the empty score cells
  for(i in 1:rows){
    for(j in 2:7){
      if(zcov3[[(i-1)*6+j-1, ncol(zcov3)]]>threshold){
        set[i, j]=1
      }
    }
  }
  set2 <- set[,-1]
  rownames(set2) <- set[,1]
  return(set)
}
makekappaset2=function(set){
  kappaset=set
  for(i in 1:nrow(set)){
    for(j in 2:7){
      if(set[[i,j]]>threshold){
        kappaset[[i,j]]=1
      }
      else{
        kappaset[[i,j]]=0
      }
    }
  }
  return(kappaset)
}
get_icc=function(set){
  subset=set[,c(2:4)]
  subset2=set[,c(5:7)]
  le=ICC(subset)
  le2=ICC(subset2)
  le=le[[1]]
  le2=le2[[1]]
  return(c(le[2,2], le2[2,2]))
}
get_kappa=function(set){
  subset=set[,c(2:4)]
  subset2=set[,c(5:7)]
  le=cohen.kappa(subset)
  le2=cohen.kappa(subset2)
  le=le[[5]]
  le2=le2[[5]]
  meankappa=(le+le2)/2
  return(meankappa)
}
#if two peaks are within a subwindow of each other, this function will remove the peak
#that is smaller in at least 5 of the 6 conditions
removepeaks=function(zcov_filtered, subwindow){
  zcov_filtered=zcov_filtered%>%mutate(remove=0)
  ncols=ncol(zcov_filtered)
  for(i in 1:nrow(zcov_filtered)){
    if(!(i==1) & zcov_filtered[[i, ncols]]==0){
      if((zcov_filtered[[i,2]]-zcov_filtered[[i-1, 2]])<=subwindow & (zcov_filtered[[i,2]]-zcov_filtered[[i-1, 2]])>0){
        if(zcov_filtered[[i-1,3]]>zcov_filtered[[i,3]]){
          zcov_filtered[[i, ncols]]=1
        }
        else{
          zcov_filtered[[i, ncols]]=-4
        }
      }
    }
    if(!(i==nrow(zcov_filtered)) & zcov_filtered[[i, ncols]]==0){
      if((zcov_filtered[[i,2]]-zcov_filtered[[i+1, 2]])>=-subwindow & (zcov_filtered[[i,2]]-zcov_filtered[[i+1, 2]])<0){
        if(zcov_filtered[[i+1,3]]>zcov_filtered[[i,3]]){
          zcov_filtered[[i, ncols]]=1
        }
        else{
          zcov_filtered[[i, ncols]]=-6
        }
      }
    }
  }
  zcov_filtered=zcov_filtered%>%group_by(position)%>%mutate(sum=sum(remove))
  zcov_filtered=zcov_filtered%>%filter(sum<=0)
  return(zcov_filtered)
}

wavelettransform=function(cov){
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  cov=wavetransform(cov)
  zcov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    mutate(score=(counts-mean(counts))/(sd(counts)))%>%mutate(originalz=(count-mean(counts))/sd(count))
  zcov_filtered=zcov%>%filter(score>9)
  zcov_filtered=removepeaks(zcov_filtered, 3)
  if(nrow(zcov_filtered)==0){
    return(0)
  }
  set=makeset(zcov, zcov_filtered)
  if(is.double(set)){
    kappaset=0
  }else{
  kappaset=makekappaset2(set)
  }
  return(list(zcov, zcov_filtered, set, kappaset))
}

#for wavelet thresholded data, scores at each nucleotide position are calculated by dividing
#the count at that region by the global count average
getwavestats=function(transcript, threshold){
  experiments=get_experiments(ribo)
  cov <- get_coverage(ribo.object = ribo,
                      name        = transcript,
                      range.lower = 29,
                      range.upper = 36,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = experiments[1:6])
  ok=wavelettransform(cov, threshold)
  if(is.double(ok)){
    return(0)
  }
  set=ok[[3]]
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  if(nrow(set)>2){
    return(get_icc(set))
  }else{
    return(0)
  }
}
getpcov=function(ribo, transcript, lengths, experiments, graph, shifts){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=lengths_aggregated(ribo)
  }
  if(missing(shifts)){
    shifts=get_pshifts(ribo, lengths, experiments, graph=FALSE)
  }
  if(missing(graph)){
    graph==FALSE
  }
  cov=get_coverage(ribo.object = ribo,
                   name        = transcript,
                   range.lower = lengths[1],
                   range.upper = lengths[1],
                   length      = FALSE,
                   alias       = TRUE,
                   compact = FALSE,
                   tidy=TRUE,
                   experiment  = experiments)
  original=cov
  for(i in 2:length(lengths)){
    temp <- get_coverage(ribo.object = ribo,
                         name        = transcript,
                         range.lower = lengths[i],
                         range.upper = lengths[i],
                         length      = FALSE,
                         alias       = TRUE,
                         compact = FALSE,
                         tidy=TRUE,
                         experiment  = experiments)
    original=rbind(original, temp)
    temp=temp%>%mutate(position=as.numeric(position))%>%arrange(desc(position))
    max=temp[[1,3]]
    temp=temp%>%mutate(position=position-shifts[[i-1]])
    temp=temp%>%arrange(position)
    
    if(shifts[[i-1]]>0){
      temp=temp%>%filter(position > 0)
      temp2=temp[c((max-shifts[[i-1]]*length(experiments)+1):max),]
      temp2=temp2%>%mutate(position=position+shifts[[i-1]])
      temp3=rbind(temp, temp2)
      cov=rbind(cov, temp3)
    }else{
      cov=rbind(cov, temp)
    }
  }
  cov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment, position)%>%summarize(count=sum(count))%>%arrange(experiment)
  original=original%>%mutate(position=as.numeric(position))%>%group_by(experiment, position)%>%summarize(count=sum(count))%>%arrange(experiment)
  if(graph==TRUE){
    original=original%>%group_by(position)%>%summarize(count=sum(count))
    pshift=cov%>%group_by(position)%>%summarize(count=sum(count))
    par(mfrow=c(1,2))
    plot(original$position, original$count, type="l")
    plot(pshift$position, pshift$count, type="l")
  }
  return(cov)
}
wavestats_p=function(transcript, threshold){
  cov=getpcov(transcript)
  ok=wavelettransform(cov, threshold)
  if(is.double(ok)){
    return(0)
  }
  set=ok[[3]]
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  if(nrow(set)>2){
    return(c(get_icc(set)))
  }else{
    return(-1)
  }
}
#clean
differentialpeaks=function(set){
  exp=ncol(set)-1
  set=set%>%mutate(statistic=0, p=0, significance=0)
  for(i in 1:nrow(set)){
    temp=data.frame("id"=c(colnames(set[,c(2:(1+exp))])), 
                    "group"=c(rep(1, exp/2), rep(2, exp/2)), 
                    "counts"=0)
    for(j in 1:exp){
      temp[j, 3]=set[i, j+1]
    }
    stats=temp%>%t_test(counts ~ group)%>%add_significance()
    set[i, 2+exp]=stats[1, 6]
    set[i, 3+exp]=stats[1, 8]
    set[i, 4+exp]=stats[1, 9]
  }
  significant=set%>%select(position, statistic, p)%>%filter(p<0.05)
  newlist=list(set, significant)
  return(newlist)
}
get_pause_sites=function(ribo, transcript, threshold, lengths, experiments){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=lengths_aggregated(ribo)
  }
  cov=get_coverage(ribo.object = ribo,
                   name        = transcript,
                   range.lower = lengths[1],
                   range.upper = lengths[length(lengths)],
                   length      = FALSE,
                   alias       = TRUE,
                   compact = FALSE,
                   tidy=TRUE,
                   experiment  = experiments)
  cov=cov%>%mutate(position=as.numeric(position))
  temp2=wavelettransform(cov, threshold)
  if(is.double(temp2)){
    return(0)
  }
  set=temp2[[3]]
  if(is.double(set)){
    return(0)
  }
  binary=temp2[[4]]
  ah=differentialpeaks(set)
  allpeaks=ah[[1]]
  allpeaks=allpeaks%>%mutate(transcript=transcript)%>%mutate(identifier=paste(position, transcript))%>%select(identifier, position, transcript, statistic, p)
  diffpeaks=ah[[2]]
  diffpeaks=diffpeaks%>%mutate(transcript=transcript)%>%mutate(identifier=paste(position, transcript))%>%select(identifier, position, transcript, statistic, p)
  binary=binary%>%mutate(transcript=transcript, identifier=paste(position, transcript))
  binary=binary[c(9, 1, 8, 2:7)]
  length=max(cov$position)
  numpeaks=nrow(allpeaks)
  diffpeaks=nrow(allpeaks%>%filter(p<0.05))
  return(list(allpeaks, binary, c(numpeaks, diffpeaks,length)))
}
get_pause_sites_l=function(ribo, transcripts, experiments, cores){
  threshold=9
  if(missing(experiments)){
  experiments=get_experiments(ribo)
  }
  lengths=lengths_aggregated(ribo)
  shifts=get_pshifts(ribo, graph=FALSE)
  pauses=mcmapply(function(a,x,y, z, w, i){
    return(get_pause_sites_l(a,x,y,z,w,i))
  }, x = transcripts, y = threshold, MoreArgs=list(a=ribo,z=shifts, w=lengths, i=experiments), mc.cores=cores)
  return(pauses)
}
#Code to run -- each function will return a 2 by 300 list with the ICC values for each
#algorithm on the 300 highest expressed transcripts in ribo. If the code throws any errors,
#please let me know and i can take a look. if it takes too long, no need to run the overlapping algo.

threshold=7
size=50
ssize=3

## CC ->  This had an error in one of the cores. This might be an edge-case. 
## It's encountered when using the top 350 transcripts. Ran fine with the top 283. 
## We might want to add some error handling to the code to isolate the problematic txn. 
transcripts=pkmtranscripts
transcripts=mmtranscripts

wavestats=mcmapply(function(x,y){
  return(getwavestats(x,y))
}, x = transcripts[1:350], y = threshold, mc.cores=3)

# zscorestats_o=mcmapply(function(x,y,z,w){
#   return(getzscorestats_o(x,y,z,w))
# }, x = transcripts, y = size, z=threshold, w=ssize, 
# mc.cores=48)


## This ran fine with top 106; Encountered_issues with 283
zscorestats_l=mcmapply(function(x,y,z,w){
  return(getzscorestats_l(x,y,z,w))
}, x = pkmtranscripts[107:108], y = size, z=threshold, w=ssize,
mc.cores=7)

#new code for the P-site adjusted coverage plots
pwavestats=mcmapply(function(x,y){
  return(wavestats_p(x,y))
}, x = transcripts[1:350], y = threshold, mc.cores=5)

listtransform=function(bad){
  length=length(bad)/2
  good=vector(mode="list", length=length)
  names=c(rep(0, length))
  for(i in 1:length){
    good[[i]]=list(bad[[i*2-1]], bad[[i*2]])
    names[i]=bad[[i*2-1]][[3]][[1]]
  }
  names(good)=names
  return(good)
}

# zscorestats_op=mcmapply(function(x,y,z,w){
#   return(zscorestats_op(x,y,z,w))
# }, x = transcripts, y = size, z=threshold, w=ssize,
# mc.cores=48)


# the peakstats functions return a list with five elements:
#1: a vector containing 14 elements, each one counting the # of pause sites detected by the algorithm in
#   a particular combination of trials. For example, the first element represents the # of pause sites that are only
#   detected in one trial, while the third element represents the # of pause sites detected in 2 trials for one condition
#   and 0 in the other. 
#2: The same as 1, but it only includes the distribution for differential pause sites.
#3: Total # of pause sites detected by the algorithm
#4: Total # of differential pause sites
#5: All nucleotide positions detected as peaks
#6: All nucleotide positions detected as differential peaks with significance and t-statistic

threshold=9
mini2=mcmapply(function(a,x,y, z, w){
  return(get_pause_sites(a,x,y,z,w))
}, x = transcripts[1:350], y = threshold, MoreArgs=list(a=ribo, z=lengths, w=experiments), mc.cores=3)

wtf=mcmapply(function(x,y,z, w, i){
  return(getpcov_l(x,y,z,w,i))
}, x = transcripts[1], w=FALSE, MoreArgs=list(lengths, experiments, shifts), mc.cores=1)


#zscorepeakcounts=mcmapply(function(x,y,z,w,a){
#  return(zscorepeakstats(x,y,z,w))
#}, x = transcripts, y = size, z=threshold, w=ssize, mc.cores=48)

# positions=set$position
# count=length(positions)
# cov_peaks=cov%>%filter(position %in% positions)
# set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
#                "control1"=1:count, "control2"=1:count, "control3"=1:count)
# set[,1]=positions
# for(i in 1:count){
#   for(j in 2:7){
#     set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)-1]
#   }
# }
# set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
#   mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))%>%mutate(sig=0)
# set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
# m=1
# for(i in 1:nrow(set)){
#   if(set[[m,8]]>5 | set[[m,9]]>5){
#     set[[m, 11]]=1
#   }
#   if(set[[m,11]]==0){
#     set=set[-m,]
#     m=m-1
#   }
#   m=m+1
# }
# set=set%>%select(-sig)
# set2 <- set[,-1]
# rownames(set2) <- set[,1]
# positions=set$position
# count=length(positions)
# pcov=cov%>%filter(position %in% positions)
# set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
#                "control1"=1:count, "control2"=1:count, "control3"=1:count)
# set[,1]=positions
# for(i in 1:count){
#   for(j in 2:7){
#     set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)]
#   }
# }
# set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
#   mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))
# set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
