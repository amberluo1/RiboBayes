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
library(numbers)
library(Metrics)

pkm.ribo <- Ribo("/Users/weiqinlu/Downloads/pkm.ribo", rename = rename_default )
mm.ribo <- Ribo("/Users/weiqinlu/Downloads/all.ribo", rename = rename_default )
all.ribo <- Ribo("/Users/amberluo/Downloads/all (1).ribo", rename = rename_default )

get_read_lengths=function(ribo, experiments, rc_threshold, ld_threshold){
  if(missing(rc_threshold)){
    rc_threshold=.85
  }
  if(missing(ld_threshold)){
    ld_threshold=.05
  }
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  rc_threshold=rc_threshold*100
  ld <- get_length_distribution(ribo.object      = ribo,
                                region      = "CDS",
                                compact=FALSE)
  num1=nrow(ld)/length(experiments)
  ld=ld%>%mutate(count=count/total.reads)
  ld=ld%>%group_by(length)%>%summarize(count=mean(count))%>%filter(count>ld_threshold)
  aggr_ld=ld$length
  rc <- get_region_counts(ribo,
                          compact=FALSE,
                          length      = FALSE,
                          transcript  = TRUE,
                          region      = c("UTR5", "UTR3","CDS"))
  rc=rc%>%group_by(experiment, length)%>%mutate(percent=count/sum(count))
  rc=rc%>%arrange(experiment)
  num=nrow(rc)/length(experiments)
  roundrc=rc%>%mutate(percent=round(percent*100, digits=0))
  aggregaterc=roundrc%>%group_by(length, region)%>%summarize(percent=mean(percent))%>%mutate(length=as.character(length))%>%mutate(experiment="all")%>%
    mutate(percent=round(percent, digits=1))
  mylengths_a=aggregaterc%>%filter(region=="CDS",percent>rc_threshold)
  aggr=mylengths_a$length
  final_aggr=intersect(aggr, aggr_ld)
  return(final_aggr)
}
plot_region_counts_by_length_e=function(ribo, lengths, experiments, expnames){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(expnames)){
    expnames=get_experiments(ribo)
  }
  if(missing(lengths)){
    ld=get_length_distribution(ribo, region="CDS", compact=FALSE)%>%filter(experiment==experiments[1])
    lengths=ld$length
  }
  rc <- get_region_counts(ribo,
                          compact=FALSE,
                          length      = FALSE,
                          transcript  = TRUE,
                          region      = c("UTR5", "UTR3","CDS"))
  rc=rc%>%filter(length %in% lengths, experiment %in% experiments)
  rc=rc%>%group_by(experiment, length)%>%mutate(percent=count/sum(count))
  rc=rc%>%arrange(experiment)
  num=nrow(rc)/length(experiments)
  shortnames=c()
  for(i in 1:length(expnames)){
    shortnames=c(shortnames, rep(expnames[i], num))
  }
  rc[,1]=shortnames
  roundrc=rc%>%mutate(percent=round(percent*100, digits=0))
  roundrc%>%ggplot(mapping=aes(x=experiment, y=percent, label=percent,fill=region))+geom_bar(position="dodge", stat="identity")+facet_wrap(~length)+
    geom_text(position=position_dodge(width=0.9), size=(12-length(expnames))*0.3, vjust=-0.25)
}
plot_region_counts_by_length=function(ribo, lengths, experiments){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    ld=get_length_distribution(ribo, region="CDS", compact=FALSE)%>%filter(experiment==experiments[1])
    lengths=ld$length
  }
  rc <- get_region_counts(ribo,
                          compact=FALSE,
                          length      = FALSE,
                          transcript  = TRUE,
                          region      = c("UTR5", "UTR3","CDS"))
  rc=rc%>%filter(length %in% lengths, experiment%in%experiments)
  rc=rc%>%group_by(experiment, length)%>%mutate(percent=count/sum(count))
  rc=rc%>%arrange(experiment)
  num=nrow(rc)/length(experiments)
  aggregaterc=rc%>%group_by(length, region)%>%summarize(percent=mean(percent))%>%mutate(length=as.character(length))%>%mutate(experiment="all")%>%
    mutate(percent=round(percent*100, digits=1))
  aggregaterc%>%ggplot(mapping=aes(x=experiment,y=percent, label=percent,fill=region))+geom_bar(position="dodge", stat="identity")+facet_wrap(~length)+
    geom_text(position=position_dodge(width=0.9), size=3, vjust=-0.25)
}
plot_metagene_by_length=function(ribo, site, positions, lengths, experiment){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=get_read_lengths(ribo)
  }
  metagene=data.frame()
  for(i in lengths){
    temp=get_metagene(ribo, site=site, range.lower=i, range.upper=i, experiment=experiment, compact=FALSE)
    temp=temp%>%mutate(length=i)
    metagene=rbind(metagene,temp)
  }
  metagene=metagene[,c(1,ncol(metagene),2:(ncol(metagene)-1))]
  colnames=colnames(metagene)
  metagene <- gather(metagene, position, counts, colnames[3]:colnames[3+radius*2], factor_key=TRUE)
  metagene=metagene%>%group_by(length)%>%mutate(counts=counts/mean(counts))
  metagene=metagene%>%group_by(length, position)%>%mutate(counts=mean(counts))%>%select(-experiment)%>%mutate(position=as.numeric(position), length=as.character(length))
  metagene%>%filter(position %in% positions)%>%ggplot(mapping=aes(x=position, y=counts, color=length))+geom_line()
}
get_pshifts=function(ribo, lengths, experiments, graph){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=get_read_lengths(ribo)
  }
  radius=metagene_radius(ribo)
  metagene=data.frame()
  for(i in lengths){
    temp=get_metagene(ribo, site="stop", range.lower=i, range.upper=i, experiment=experiments, compact=FALSE)
    temp=temp%>%mutate(length=i)
    metagene=rbind(metagene,temp)
  }
  metagene=metagene[,c(1,ncol(metagene),2:(ncol(metagene)-1))]
  colnames=colnames(metagene)
  metagene <- gather(metagene, position, counts, colnames[3]:colnames[3+radius*2], factor_key=TRUE)
  metagene=metagene%>%group_by(length)%>%mutate(counts=counts/mean(counts))
  metagene=metagene%>%group_by(length, position)%>%mutate(counts=mean(counts))%>%filter(experiment %in% experiments[1])%>%select(-experiment)%>%mutate(position=as.numeric(position))
  mod0=c(seq(3, 3+radius, by=3))
  mod1=c(seq(1, 1+radius, by=3))
  mod2=c(seq(2, 2+radius, by=3))
  metagene0=metagene%>%filter(position %in% mod0)%>%group_by(length)%>%summarize(mean=mean(counts))%>%mutate(mod=0)
  metagene1=metagene%>%filter(position %in% mod1)%>%group_by(length)%>%summarize(mean=mean(counts))%>%mutate(mod=1)
  metagene2=metagene%>%filter(position %in% mod2)%>%group_by(length)%>%summarize(mean=mean(counts))%>%mutate(mod=2)
  mods=rbind(metagene0, metagene1, metagene2)
  data_wide=spread(mods, mod, mean)
  data_wide = data_wide%>%mutate(max=ifelse(`0`>`1` & `0`>`2`, 0, ifelse(`1`>`2`, 1, 2)))
  base=data_wide[[1,5]]
  shifts=c()
  a=metagene%>%filter(length %in% lengths[1])
  pcov=a
  for(i in 2:nrow(data_wide)){
    min=500000
    max=data_wide[[i,data_wide[[i,5]]+2]]
    for(inc in 0:2){
      b=metagene%>%filter(length %in% lengths[i])
      bcounts=b[c(1:(nrow(b)-inc)),]$counts
      bcounts=c(rep(0, inc), bcounts)
      #acounts=a[c((inc+1):nrow(a)),]$counts
      mod=(base-inc) %% 3
      if(rmse(a$counts, bcounts)*max/data_wide[[i,mod+2]] < min){
        min=rmse(a$counts, bcounts)*max/data_wide[[i,mod+2]]
        shift=inc
      }
    }
    bcounts=b[c(1:(nrow(b)-shift)),]$counts
    bcounts=c(rep(0, shift), bcounts)
    pcov=rbind(pcov, data.frame("length"=rep(lengths[i], nrow(b)), "position"=c(1:nrow(b)), "counts"=bcounts))
    shifts=c(shifts, shift)
  }
  myacf=pcov%>%group_by(position)%>%summarize(counts=mean(counts))
  originalacf=metagene%>%group_by(position)%>%summarize(counts=mean(counts))
  names(shifts)=lengths[c(2:length(lengths))]
  if(graph==TRUE){
    original=plot_metagene_by_length(ribo, "stop", positions=c(1:get_info(ribo)$attributes$left_span), lengths=lengths, experiments=experiments)
    pshift=pcov%>%mutate(length=as.character(length))%>%filter(position %in% c(1:get_info(ribo)$attributes$left_span))%>%ggplot(mapping=aes(x=position, y=counts, color=length))+geom_line()
    print(original + pshift)
    metagene=data.frame()
    for(i in lengths){
      temp=get_metagene(ribo, site="stop", range.lower=i, range.upper=i, experiment=experiments, compact=FALSE)
      temp=temp%>%mutate(length=i)
      metagene=rbind(metagene,temp)
    }
    metagene=metagene[,c(1,ncol(metagene),2:(ncol(metagene)-1))]
    colnames=colnames(metagene)
    metagene <- gather(metagene, position, counts, colnames[3]:colnames[3+radius*2], factor_key=TRUE)
    metagene=metagene%>%group_by(length)%>%mutate(counts=counts/mean(counts))
    metagene=metagene%>%group_by(position)%>%mutate(counts=sum(counts)/length(experiments))%>%mutate(position=as.numeric(position))
    original=metagene%>%filter(position %in% c(1:get_info(ribo)$attributes$left_span))%>%ggplot(mapping=aes(x=position, y=counts))+geom_line()
    pshift=pcov%>%group_by(position)%>%mutate(counts=sum(counts))%>%filter(position %in% c(1:get_info(ribo)$attributes$left_span))%>%ggplot(mapping=aes(x=position, y=counts))+geom_line()
    print(original+ pshift)
  }
  return(shifts)
}
plot_coverage_by_length=function(ribo, transcript, lengths, positions){
  cov=get_coverage(ribo.object = ribo,
                   name        = transcript,
                   range.lower = lengths[1],
                   range.upper = lengths[length(lengths)],
                   length      = FALSE,
                   alias       = TRUE,
                   compact = FALSE,
                   tidy=TRUE,
                   experiment  = experiments)
  cov=cov%>%mutate(position=as.numeric(position), length=as.character(length))
  cov=cov%>%group_by(position, length)%>%summarize(count=sum(count))
  cov=cov%>%filter(position %in% positions)
  ah=cov%>%ggplot(mapping=aes(x=position, y=count, color=length))+geom_line()
  print(ah)
}

#clean - may revise
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
  if(nrow(set)==0){
    return(0)
  }
  numexperiments=ncol(set)-1
  set=set%>%mutate(avg1=rowMeans(set[,c(2:(numexperiments/2+1))]), avg2=rowMeans(set[,c((numexperiments/2+2):(numexperiments+1))]))
  m=1
  for(i in 1:nrow(set)){
    if(set[[m,(2+numexperiments)]]<6 & set[[m,(3+numexperiments)]]<6){
      set=set[-m,]
      m=m-1
    }
    m=m+1
  }
  if(nrow(set)==0){
    return(0)
  }
  set=set%>%select(-avg1, -avg2)
  positions=set$position
  mypositions=zcov%>%group_by(experiment)%>%mutate(originalz=(count-mean(count))/sd(count))%>%filter(position %in% positions)
  originalset=dcast(mypositions, position ~ experiment, value.var="originalz")
  return(list(set, originalset))
}
makekappaset=function(zcov, zcov_filtered){
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
#clean
makekappaset2=function(set){
  kappaset=set
  for(i in 1:nrow(set)){
    for(j in 2:7){
      if(set[[i,j]]>6){
        kappaset[[i,j]]=1
      }
      else{
        kappaset[[i,j]]=0
      }
    }
  }
  return(kappaset)
}
#clean
get_icc=function(set){
  subset=set[,c(2:4)]
  subset2=set[,c(5:7)]
  le=ICC(subset)
  le2=ICC(subset2)
  le=le[[1]]
  le2=le2[[1]]
  return(c(le[2,2], le2[2,2]))
}
#clean
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
getpcov=function(ribo, transcript, lengths, experiments, graph, shifts){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=get_read_lengths(ribo)
  }
  if(missing(shifts)){
    shifts=get_pshifts(ribo, lengths, experiments, graph=FALSE)
  }
  if(missing(graph)){
    graph==TRUE
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
    temp=temp%>%mutate(position=as.numeric(position)+shifts[[i-1]])
    temp=temp%>%arrange(position)
    if(shifts[[i-1]]>0){
      temp=temp%>%filter(position <= max)
      temp2=temp[c(1:(shifts[[i-1]]*length(experiments))),]
      temp2=temp2%>%mutate(position=position-shifts[[i-1]])
      temp3=rbind(temp2, temp)
      cov=rbind(temp3, cov)
    }else{
      cov=rbind(cov, temp)
    }
  }
  ncov=cov%>%mutate(position=as.numeric(position))
  cov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment, position)%>%summarize(count=sum(count))%>%arrange(experiment)
  original=original%>%mutate(position=as.numeric(position))%>%group_by(experiment, position)%>%summarize(count=sum(count))%>%arrange(experiment)
  if(graph==TRUE){
    original=original%>%group_by(position)%>%summarize(count=sum(count))
    pshift=cov%>%group_by(position)%>%summarize(count=sum(count))
    par(mfrow=c(1,2))
    plot(original$position, original$count, type="l")
    plot(pshift$position, pshift$count, type="l")
    par(mfrow=c(1,1))
    transcript2=transcript
    rc=get_internal_region_coordinates(pkm.ribo, alias=TRUE)%>%select(UTR3_start, transcript)%>%filter(transcript==transcript2)
    stop=rc[[1,1]]
    ocov=get_coverage(ribo.object = ribo,
                      name        = transcript,
                      range.lower = lengths[1],
                      range.upper = lengths[length(lengths)],
                      length      = FALSE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = experiments)
    ocov=ocov%>%mutate(position=as.numeric(position), length=as.character(length))
    ocov=ocov%>%group_by(position, length)%>%summarize(count=sum(count))
    original=ocov%>%filter(position %in% c((stop-120):(stop-80)))%>%ggplot(mapping=aes(x=position, y=count, color=length))+geom_line()
    pshift=ncov%>%group_by(position, length)%>%summarize(count=sum(count))
    pshift=pshift%>%mutate(position=as.numeric(position), length=as.character(length))%>%filter(position %in% c((stop-120):(stop-80)))
    pshift=pshift%>%ggplot(mapping=aes(x=position, y=count, color=length))+geom_line()
    print(original + pshift)
    wow=ocov%>%group_by(position)%>%summarize(count=sum(count))
    wow2=ncov%>%group_by(position)%>%summarize(count=sum(count))
    ocov=ocov%>%group_by(position)%>%summarize(count=sum(count))%>%mutate(position=as.numeric(position))
    ncov=ncov%>%group_by(position)%>%summarize(count=sum(count))%>%mutate(position=as.numeric(position))
    ocov2=ocov%>%filter(position %in% c((stop-50):(stop-20)))
    ncov2=ncov%>%filter(position %in% c((stop-50):(stop-20)))
    original=ocov2%>%ggplot(mapping=aes(x=position, y=count))+geom_line()
    pshift=ncov2%>%ggplot(mapping=aes(x=position, y=count))+geom_line()
    print(original + pshift)
  }
  return(cov)
}
#clean
wavetransform=function(cov){
  experiments=unique(cov$experiment)
  max=max(cov$position)
  power=floor(log(max, 2))
  noise=c(rep(0.001,nrow(cov)))
  noise=jitter(noise, factor = 1, amount=NULL)
  wavecov=cov%>%arrange(experiment)%>%ungroup()%>%mutate(noises=noise, ocounts=counts)
  transformed=c()
  for(i in 1:length(experiments)){
    subcov=wavecov%>%filter(experiment %in% experiments[i])%>%mutate(counts=counts+noises)
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
  wavecov=wavecov[,-5]
  return(wavecov)
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
adjustscore=function(zcov){
  zcov=zcov%>%arrange(position)
  max=max(zcov$position)
  numexperiments=nrow(zcov)/max
  zcov=zcov%>%mutate(localmean=0, localsd=0, group=0)
  ceiling=ceiling(max/101)
  allmeans=c()
  allsds=c()
  singlemeans=c()
  singlesds=c()
  for(i in 1:2){
    for(j in 1:ceiling){
      zcov1=zcov[(c((i-1)*nrow(zcov)/2)+1):(nrow(zcov)/(3-i)),]
      subcov=zcov1%>%filter(ceiling(position/101)==j)
      mean=mean(subcov$counts)
      sd=sd(subcov$counts)
      singlemeans=c(singlemeans, mean)
      singlesds=c(singlesds, sd)
      subcov=subcov%>%mutate(localmean=mean, localsd=sd)
      localmean=subcov$localmean
      localsd=subcov$localsd
      allmeans=c(allmeans, localmean)
      allsds=c(allsds, localsd)
    }
  }
  cov=zcov
  cov[,8]=allmeans
  cov[,9]=allsds
  cov=cov%>%mutate(group=ceiling(position/101),adjmean=0, adjsd=0)
  cov=cov%>%group_by(group)%>%mutate(percent=(position-group*101+50)/(100))
  cov2=cov%>%filter(group>1 & group < length(singlemeans)/2)
  cov2=cov2%>%filter(score>4)%>%mutate(newscore=ifelse(percent==0, score, ifelse(percent<0, abs(percent)*(counts-singlemeans[group-1])/(singlesds[group-1])+((1-abs(percent))*(counts-localmean)/localsd),
                                                                                 ifelse(percent>0, percent*(counts-singlemeans[group+1])/(singlesds[group+1])+((1-percent)*(counts-localmean)/localsd), 0))))
  cov=cov%>%filter(score<=4)
  cov2=cov2%>%mutate(score=newscore)%>%select(-newscore)
  cov=rbind(cov, cov2)
  cov=cov%>%arrange(position)
  return(cov)
}
#clean
wavelettransform=function(cov){
  no0=cov%>%filter(count>0)
  no0=no0%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  cov=cov%>%filter(count==0)%>%mutate(counts=0)
  cov=rbind(cov, no0)
  cov=cov%>%arrange(position)
  no0=cov%>%filter(count>0)
  osd=sd(no0$counts)
  omean=mean(no0$counts)
  cov=wavetransform(cov)
  no0=cov%>%filter(count>0)
  mean=mean(no0$counts)
  sd=sd(no0$counts)
  zcov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    mutate(score=(counts-mean)/(sd))%>%mutate(originalz=(ocounts-omean)/osd)
  zcov_filtered=zcov%>%filter(score>7)
  if(nrow(zcov_filtered)==0){
    return(0)
  }
  zcov_filtered=removepeaks(zcov_filtered, 3)
  if(nrow(zcov_filtered)==0){
    return(0)
  }
  sets=makeset(zcov, zcov_filtered)
  if(is.double(sets)){
    return(0)
  }
  set=sets[[1]]
  originalset=sets[[2]]
  if(is.double(set)){
    kappaset=0
  }else{
    kappaset=makekappaset2(set)
  }
  return(list(zcov, zcov_filtered, originalset, kappaset))
}
#clean
get_pause_sites=function(ribo, transcript, shifts, lengths, experiments){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=get_read_lengths(ribo)
  }
  cov=getpcov(ribo, transcript, lengths, experiments, graph=FALSE, shifts)
  cov=cov%>%mutate(position=as.numeric(position))
  temp=wavelettransform(cov)
  if(is.double(temp)){
    return(0)
  }
  originalset=temp[[3]]
  if(is.double(originalset)){
    return(0)
  }
  binary=temp[[4]]
  ah=differentialpeaks(originalset)
  allpeaks=ah[[1]]
  allpeaks=allpeaks%>%mutate(transcript=transcript)%>%mutate(identifier=paste(position, transcript))%>%select(identifier, position, transcript, statistic, p)
  diffpeaks=ah[[2]]
  diffpeaks=diffpeaks%>%mutate(transcript=transcript)%>%mutate(identifier=paste(position, transcript))%>%select(identifier, position, transcript, statistic, p)
  binary=binary%>%mutate(transcript=transcript, identifier=paste(position, transcript))
  originalset=originalset%>%mutate(transcript=transcript, identifier=paste(position, transcript))
  binary=binary[c(9, 1, 8, 2:7)]
  originalset=originalset[c(9, 1, 8, 2:7)]
  length=max(cov$position)
  numpeaks=nrow(allpeaks)
  diffpeaks=nrow(allpeaks%>%filter(p<0.05))
  return(list(allpeaks, originalset, binary, c(numpeaks, diffpeaks,length)))
}
#temporary function for bad bad situations :(
listtransform=function(bad){
  length=length(bad)/3
  good=vector(mode="list", length=length)
  names=c(rep(0, length))
  for(i in 1:length){
    good[[i]]=list(bad[[i*3-2]], bad[[i*3-1]], bad[[i*3]])
    names[i]=bad[[i*3-2]][[3]][[1]]
  }
  names(good)=names
  return(good)
}

##################################################
ribo=pkm.ribo
lengths=get_read_lengths(ribo)
experiments=get_experiments(ribo)
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
high_exp = rowSums( cpm(rc_CDS_w[,-1]) > 369) > 1
sum(high_exp)

transcripts=rc_CDS_w$transcript[high_exp]

sites=mcmapply(function(a,x,y, z, w){
  return(get_pause_sites(a,x,y,z,w))
}, x = transcripts, MoreArgs=list(a=ribo, y=shifts, z=lengths, w=experiments), mc.cores=3)
