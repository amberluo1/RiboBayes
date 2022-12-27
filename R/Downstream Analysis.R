library(ggvenn)
library(UpSetR)
removezeros=function(sites){
  zeropos=c()
  m=1
  n=1
  for(i in 1:length(sites)){
    if(is.double(sites[[n]])){
      zeropos[m]=names(sites[n])
      sites=sites[-n]
      m=m+1
      n=n-1
    }
    n=n+1
  }
  return(sites)
}
get_ps_info=function(sites){
  sites=removezeros(sites)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[2]])
  }
  return(data)
}
get_ps_change=function(sites, bayes.adjust){
  if(missing(bayes.adjust)){
    bayes.adjust=FALSE
  }
  sites=removezeros(sites)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[1]])
  }
  return(data)
}
plot_pause_site_conservation=function(sites){
  sites=removezeros(sites)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[3]])
  }
  namedata=data
  for(i in 4:(4+length(experiments)-1)){
    for(j in 1:nrow(data)){
      if(data[[j,i]]==1){
        namedata[j,i]=data[[j,1]]
      }else{
        namedata[j,i]=0
      }
    }
  }
  numexperiments=(ncol(namedata)-3)/2
  categories=list()
  m=1
  for(i in 1:2){
    for(j in 1:numexperiments){
      scripts=as.list(namedata[,m+3])
      n=1
      for(k in 1:length(scripts)){
        if(scripts[[n]]==0){
          scripts=scripts[-n]
          n=n-1
        }
        n=n+1
      }
      categories[[m]]=as.list(scripts)
      m=m+1
    }
  }
  experiments=colnames(data)[4:ncol(data)]
  names(categories)=experiments
  exp1=categories[1:numexperiments]
  exp2=categories[(numexperiments+1):(numexperiments*2)]
  first=ggvenn(
    exp1,
    fill_color = c("45534CAA", "#3BC000FF", "#0096FF"),
    stroke_size = 1, set_name_size = 4
  )
  second=ggvenn(
    exp2,
    fill_color = c("45534CAA", "#3BC000FF", "#0096FF"),
    stroke_size = 1, set_name_size = 4
  )
  print(first+second)
  data=data%>%mutate_at(vars(experiments), funs(as.numeric(.)))
  print(upset(data, nsets = 6, nintersects = NA, mb.ratio = c(0.5, 0.5),keep.order=TRUE,sets=colnames(data)[4:9],
              order.by = "freq", decreasing = c(TRUE,FALSE)))
}
plot_pause_site_regulation=function(sites){
  sites=removezeros(sites)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[1]])
  }
  data=data%>%mutate("Pause Site Expression"=ifelse(p<0.05, "Differential", "Constant"), statistic=statistic * -1)
  print(ggstatsplot::ggbetweenstats(
    data = data,
    x = `Pause Site Expression`,
    y = statistic,
    outlier.tagging = TRUE, outlier.coef=2.0,
    outlier.label = transcript,
    title = "Distribution of Pause Sites in 350 Highly Expressed Transcripts from eIF5A ovx and PKM KD Cells",
    ylab = "t-statistic of Pause Site Expression Change in eIF5A ovx vs. control "
  )
  )
}
plot_pause_site_distribution=function(ribo, sites){
  sites=removezeros(sites)
  rc=get_internal_region_coordinates(ribo, alias=TRUE)
  lengths=get_internal_region_lengths(ribo, alias=TRUE)
  rc=rc%>%left_join(lengths)%>%mutate(startcodon=CDS_start-10, length=CDS)%>%select(transcript, startcodon, length)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[1]])
  }
  data=data%>%left_join(rc)
  data=data%>%mutate(relativepos=(position-startcodon)/length, distance=position-startcodon)%>%filter(distance>0 & relativepos<1)
  data=data%>%mutate(sig=ifelse(p<0.05, "Differential", "Constant"))
  sigdata=data%>%filter(p<0.05)
  a=data%>%ggplot(mapping=aes(x=relativepos, fill=sig))+geom_density(alpha=0.4)
  sigdata=sigdata%>%mutate(regulation=ifelse(statistic>0, "Downregulated", "Upregulated"))
  b=sigdata%>%ggplot(mapping=aes(x=relativepos, fill=regulation))+geom_density(alpha=0.4)
  print(a+b)
}
plot_pause_sites_by_transcript=function(experiment_name, sites){
  if(missing(experiment_name)){
    experiment_name="ribo"
  }
  sites=removezeros(sites)
  length=length(sites)
  numpauses=data.frame("transcript"=names(sites), "ps"=1:length, "dps"=1:length, "length"=1:length)

  ps=mcmapply(function(i){
    return(sites[[i]][[4]][[1]])
  }, i=1:length(sites), mc.cores=3)

  dps=mcmapply(function(i){
    return(sites[[i]][[4]][[2]])
  }, i=1:length(sites), mc.cores=3)

  length=mcmapply(function(i){
    return(sites[[i]][[4]][[3]])
  }, i=1:length(sites), mc.cores=3)

  numpauses[2]=ps
  numpauses[3]=dps
  numpauses[4]=length
  numpauses=numpauses%>%mutate(ribo=experiment_name)

  print(ggstatsplot::ggbetweenstats(
    data = numpauses,
    x = ribo,
    y = ps,
    outlier.tagging = TRUE, outlier.coef=2.0,
    outlier.label = transcript,
    title = "Number of Pause Sites by Transcript",
    xlab = experiment_name,
    ylab = "Number of Constant & Differential Pause Sites"
  ))
  dps_t=table(dps)
  par(mfrow=c(1,2))
  bp1=barplot(dps_t, col=c("lightblue"), main="Number of Differential Pause Sites per Transcript", xlab="Number of Differential Pause Sites", ylab="Number of Transcripts")
  text(bp1,0, round(dps_t, 1),cex=0.9,pos=3)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[1]])
  }
  updata=data%>%filter(p<0.05 & statistic<0)
  up=plyr::count(updata, "transcript")
  up=up%>%dplyr::rename(up=freq)
  downdata=data%>%filter(p<0.05 & statistic>0)
  down=plyr::count(downdata, "transcript")
  down=down%>%dplyr::rename(down=freq)
  zerodata=data.frame("transcript"=names(sites), "up"=rep(0, length(sites)), "down"=rep(0, length(sites)))
  zerodata=zerodata%>%rows_update(up)
  zerodata=zerodata%>%rows_update(down)
  Upregulated=table(zerodata$up)
  Downregulated=table(zerodata$down)
  full=rbind(Upregulated, Downregulated)
  if(length(Upregulated)!=length(Downregulated)){
    if(length(Upregulated)<length(Downregulated)){
      diff=length(Downregulated)-length(Upregulated)
      for(i in 1:diff){
        full[[length(Upregulated)*2+1+(i-1)*2]]=0
      }
    }
  }
  y=full[1:length(full)]
  max=max(y)
  bp2=barplot(full, col=c("lightblue", "pink"), beside=TRUE,
              main="Upregulated & Downregulated Pause Sites by Transcript", xlab="Differential Pause Sites by Expression Change", ylab="Number of Transcripts", legend=c("Upregulated", "Downregulated"),
              args.legend=list(title="Pause Sites", pt.cex=1.3,cex=0.7,"topright", text.font=3),ylim=c(0,ceiling(max/100)*100))
  text(bp2,y+7,labels=as.character(y), cex=0.8)
}
