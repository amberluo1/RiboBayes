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
get_pshifts=function(ribo, lengths, experiments, graph){
  if(missing(experiments)){
    experiments=get_experiments(ribo)
  }
  if(missing(lengths)){
    lengths=get_read_lengths(ribo)
  }
  radius=metagene_radius(ribo)
  metagene=data.frame()
  for(i in as.numeric(lengths)){
    temp=get_metagene(ribo, site="stop", range.lower=i, range.upper=i, experiment=experiments, compact=FALSE)
    temp=temp%>%mutate(length=i)
    metagene=rbind(metagene,temp)
  }
  metagene=metagene[,c(1,ncol(metagene),2:(ncol(metagene)-1))]
  colnames=colnames(metagene)
  metagene <- gather(metagene, position, counts, colnames[3]:colnames[3+radius*2], factor_key=TRUE)
  metagene=metagene%>%group_by(length)%>%mutate(counts=counts/mean(counts))
  metagene=metagene%>%group_by(length, position)%>%mutate(counts=mean(counts))%>%filter(experiment %in% experiments[1])%>%dplyr::select(-experiment)%>%mutate(position=as.numeric(position))
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
    pcov=rbind(pcov, data.frame("length"=rep(as.numeric(lengths)[i], nrow(b)), "position"=c(1:nrow(b)), "counts"=bcounts))
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
    rc=get_internal_region_coordinates(pkm.ribo, alias=TRUE)%>%dplyr::select(UTR3_start, transcript)%>%filter(transcript==transcript2)
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
