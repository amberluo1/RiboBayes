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
  significant=set%>%dplyr::select(position, statistic, p)%>%filter(p<0.05)
  newlist=list(set, significant)
  return(newlist)
}
bayesian_p_adjust=function(sites, graph){
  if(missing(graph)){
    graph=FALSE
  }
  numexperiments=(ncol(sites[[1]][[2]])-3)/2
  exptype=c()
  vector=c("CTR", "EXP")
  for(i in 1:2){
    for(j in 1:numexperiments){
      exptype=c(exptype, vector[i])
    }
  }
  exptype <- factor(exptype)
  counts=get_ps_info(sites)
  counts=counts%>%dplyr::select(-position, -transcript)
  for(i in 2:ncol(counts)){
    counts=counts[counts[,i]>=0,]
  }
  identifiers=counts$identifier
  counts2=get_ps_info(sites)
  counts2=counts2%>%filter(!(identifier%in%(identifiers)))%>%dplyr::select(-position, -transcript)
  counts2[counts2 < 0] <- 0
  counts=rbind(counts, counts2)
  y <- DGEList(counts=counts[,-1],group=exptype, genes = counts[,1])
  y <- calcNormFactors(y, method = "TMM")

  design <- model.matrix(~0 + exptype)
  colnames(design) <- levels(exptype)
  y <- estimateDisp(y,design)
  if(graph==TRUE){
    plotBCV(y)
    plotMDS(y)
  }
  fit <- glmQLFit(y,design)
  my.contrasts <- makeContrasts(
    PS_CHANGE = EXP-CTR,
    levels = design
  )
  qlfChange <- glmQLFTest(fit, contrast=my.contrasts[,"PS_CHANGE"])
  qlfChange=as.data.frame(qlfChange)
  qlfChange=qlfChange%>%rename(identifier=genes, p=PValue)%>%dplyr::select(identifier, p, logFC)%>%
    mutate(position=as.numeric(sub(" .*", "", identifier)), transcript=sub(".* ", "", identifier))
  if(length(get_experiments(ribo)>6)){
    qlfChange=qlfChange[,c(1, 4, 5, 2, 3)]%>%rename(statistic=logFC)%>%mutate(statistic=-statistic)
  }
  return(qlfChange)
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
  cov2=cov2%>%mutate(score=newscore)%>%dplyr::select(-newscore)
  cov=rbind(cov, cov2)
  cov=cov%>%arrange(position)
  return(cov)
}
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
get_pause_sites_h=function(ribo, transcript, shifts, lengths, experiments, cores){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=get_read_lengths(ribo)
  }
  if(missing(shifts)){
    shifts=get_pshifts(ribo, lengths, experiments, FALSE)
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
  allpeaks=allpeaks%>%mutate(transcript=transcript)%>%mutate(identifier=paste(position, transcript))%>%dplyr::select(identifier, position, transcript, statistic, p)
  diffpeaks=ah[[2]]
  diffpeaks=diffpeaks%>%mutate(transcript=transcript)%>%mutate(identifier=paste(position, transcript))%>%dplyr::select(identifier, position, transcript, statistic, p)
  binary=binary%>%mutate(transcript=transcript, identifier=paste(position, transcript))
  originalset=originalset%>%mutate(transcript=transcript, identifier=paste(position, transcript))
  binary=binary[c(9, 1, 8, 2:7)]
  originalset=originalset[c(9, 1, 8, 2:7)]
  length=max(cov$position)
  numpeaks=nrow(allpeaks)
  diffpeaks=nrow(allpeaks%>%filter(p<0.05))

  return(list(allpeaks, originalset, binary, c(numpeaks, diffpeaks,length)))
}

#' Create a RiboBayes object from a Ribo object
#'
#' Takes as input a Ribo object and returns a RiboBayes object containing information on the
#' location and differential expression data of all ribosome pause sites on the inputted transcripts.
#'
#' @param ribo Ribo object
#' @param transcripts Character vector of transcript names.
#' @param p_shift if `FALSE`, P-site adjustments are not applied. Defaults to `TRUE`.
#' @param experiments Character vector of ribosome profiling experiments.
#' @param cores Number of cores to run `get_pause_sites()` on.
#' @return A RiboBayes object containing []
#' @export
get_pause_sites=function(ribo, transcripts, p_shift, experiments, cores){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(cores)){
    cores=1
  }
  if(missing(p_shift)){
    p_shift=TRUE
  }
  if(p_shift){
    lengths=get_read_lengths(ribo)
    shifts=get_pshifts(ribo, lengths, experiments, graph=FALSE)
  }else{
    lengths=(length_min(ribo):length_max(ribo))
    temp=length_max(ribo)-length_min(ribo)+1
    shifts=integer(temp)
    names(shifts)=lengths
  }
  sites=mcmapply(function(a,x,y, z, w){
    return(get_pause_sites_h(a,x,y,z,w))
  }, x = transcripts, MoreArgs=list(a=ribo, y=shifts, z=lengths, w=experiments), mc.cores=cores)
  return(sites)
}
