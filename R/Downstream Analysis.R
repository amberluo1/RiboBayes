library(ggvenn)
pausedata=function(sites){
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[1]])
  }
  return(data)
}
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
}t
transcriptdata=function(sites){
  data=data.frame("transcript"=1:length(sites), "numpeaks"=1:length(sites), "length"=1:length(sites))
  data[,1]=names(sites)
  for(i in 1:length(sites)){
  data[i,2]=nrow(sites[[i]][[1]])
  data[i,3]=sites[[i]][[3]]
  }
  return(data)
}
plot_pause_site_conservation=function(sites){
  sites=removezeros(sites)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[2]])
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
oops=function(data){
  for(i in 4:9){
    for(j in 1:nrow(data)){
      if(data[[j, i]]!=0){
        data[[j, i]]=1
      }
    }
  }
  return(data)
}
plot_pause_site_regulation=function(sites){
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[1]])
  }
  data=data%>%mutate("Pause Site Expression"=ifelse(p<0.05, "Differential", "Constant"))
  print(ggstatsplot::ggbetweenstats(
    data = data,
    x = "Pause Site Expression",
    y = statistic,
    outlier.tagging = TRUE, outlier.coef=1.5,
    outlier.label = transcript,
    title = "Distribution of Pause Sites in 350 Highly Expressed Transcripts from eIF5A ovx and PKM KD Cells",
    ylab = "Number of Pause Sites by Transcript"
  )
)
}
