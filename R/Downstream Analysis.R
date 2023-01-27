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

#' Create a data frame summarizing the output of `get_pause_sites()`
#'
#' Takes as input a list of ribosome pause sites (the output of `get_pause_sites()`) and returns a
#' data frame with one row for each pause site, with columns containing an unique identifier, transcript location,
#' nucleotide position, and significance score in all replicates for each pause site.
#'
#' @param sites List of ribosome pause sites; usually the output of `get_pause_sites()`
#' @return Data frame with one row for each pause site and columns containing a unique identifier,
#' transcript location, nucleotide position, and significance score in all replicates for each pause site.
#' @export
get_ps_info=function(sites){
  sites=removezeros(sites)
  data=data.frame()
  for(i in 1:length(sites)){
    data=rbind(data, sites[[i]][[2]])
  }
  return(data)
}

#' Create a data frame with the probability of differential expression for each pause site
#'
#' Takes as input a list of ribosome pause sites (the output of `get_pause_sites()`) and returns a
#' data frame with one row for each pause site and columns containing information about pause site location
#' and differential expression (p-value, logFC).
#'
#' @param sites List of ribosome pause sites; usually the output of `get_pause_sites()`
#' @return Data frame with one row for each pause site and columns containing a unique identifier,
#' transcript location, nucleotide position, log fold change (logFC), and p-value of change across conditions for each pause site.
#' @details The function \code{\link{get_ps_change}} calls edgeR's methods to quantify change across conditions
#' with a limited number of samples, producing a log fold-change (logFC) and p-value. You can read more about these methods
#' \href{https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf}{here}.
#' @export
get_ps_change=function(sites){
  sites=removezeros(sites)
  # data=data.frame()
  # for(i in 1:length(sites)){
  #   data=rbind(data, sites[[i]][[1]])
  # }
  data=bayesian_p_adjust(sites)
  return(data)
}

#' Create Venn Diagrams to visualize the conservation of pause sites across experiments
#'
#' Takes as input a list of ribosome pause sites (the output of `get_pause_sites()`) and plots
#' Venn diagrams depicting the overlap of detected pause sites across experiments.
#'
#' @param sites List of ribosome pause sites; usually the output of `get_pause_sites()`
#' @return One figure with two Venn diagrams (one for each condition) depicting the overlap of pause site
#' expression across experiments. One UpSetR graph providing more detailed information about the expression
#' of pause sites across replicates and experiments.
#' @export
plot_ps_conservation=function(sites){
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

#' Plot violin plots to visualize the distribution of t-statistics for pause site significance.
#'
#' Takes as input a list of ribosome pause sites (the output of `get_pause_sites()`) and plots
#' two violin plots showing the distribution of t-statistics for constant and differential pause sites.
#'
#' @param sites List of ribosome pause sites; usually the output of `get_pause_sites()`
#' @return A `ggstatsplot` figure with two violin plots: one showing the distribution of pause site
#' t-statistics for constant pause sites only (p>0.05) and one showing the distribution of pause
#' site t-statistics of differential pause sites only (p<=0.05).
#' @export
plot_ps_regulation=function(sites){
  sites=removezeros(sites)
  data=bayesian_p_adjust(sites)
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

#' Create Venn Diagrams to visualize the conservation of pause sites across experiments
#'
#' Takes as input a list of ribosome pause sites (the output of `get_pause_sites()`) and plots
#' Venn diagrams depicting the overlap of detected pause sites across experiments.
#'
#' @param sites List of ribosome pause sites; usually the output of `get_pause_sites()`
#' @return One figure with two Venn diagrams (one for each condition) depicting the overlap of pause site
#' expression across experiments. One UpSetR graph providing more detailed information about the expression
#' of pause sites across replicates and experiments.
#' @export
plot_ps_distribution=function(ribo, sites){
  sites=removezeros(sites)
  rc=get_internal_region_coordinates(ribo, alias=TRUE)
  lengths=get_internal_region_lengths(ribo, alias=TRUE)
  rc=rc%>%left_join(lengths)%>%mutate(startcodon=CDS_start-10, length=CDS)%>%select(transcript, startcodon, length)
  data=bayesian_p_adjust(sites)
  data=data%>%left_join(rc)
  data=data%>%mutate(relativepos=(position-startcodon)/length, distance=position-startcodon)%>%filter(distance>0 & relativepos<1)
  data=data%>%mutate(sig=ifelse(p<0.05, "Differential", "Constant"))
  sigdata=data%>%filter(p<0.05)
  a=data%>%ggplot(mapping=aes(x=relativepos, fill=sig))+geom_density(alpha=0.4)
  sigdata=sigdata%>%mutate(regulation=ifelse(statistic>0, "Downregulated", "Upregulated"))
  b=sigdata%>%ggplot(mapping=aes(x=relativepos, fill=regulation))+geom_density(alpha=0.4)
  print(a+b)
}

#' Plot bar plots to visualize the distribution of differential pause sites per transcript
#'
#' Takes as input a list of ribosome pause sites (the output of `get_pause_sites()`) and plots
#' two bar plots: one showing the distribution of the number of differential pause sites per transcript,
#' and another showing the distribution of the number of upregulated vs. downregulated pause sites per
#' transcript.
#'
#' @param sites List of ribosome pause sites; usually the output of `get_pause_sites()`
#' @return `NULL`
#' @export
plot_dps_per_transcript=function(sites){
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

  dps_t=table(dps)
  par(mfrow=c(1,2))
  bp1=barplot(dps_t, col=c("gray"), main="Number of Differential Pause Sites per Transcript", xlab="Number of Differential Pause Sites", ylab="Number of Transcripts")
  text(bp1,0, round(dps_t, 1),cex=0.9,pos=3)
  data=get_ps_change(sites)

  # extract upregulated and downregulated pause sites
  # tally up number of differential pause sites by transcript

  updata=data%>%filter(p<0.05 & statistic<0)
  up=plyr::count(updata, "transcript")%>%dplyr::rename(up=freq)

  downdata=data%>%filter(p<0.05 & statistic>0)
  down=plyr::count(downdata, "transcript")%>%dplyr::rename(down=freq)

  # create data frame with all transcripts, regardless of whether they have DPS or not

  all_transcripts=data.frame("transcript"=names(sites), "up"=rep(0, length(sites)), "down"=rep(0, length(sites)))
  all_transcripts=all_transcripts%>%rows_update(up)
  all_transcripts=all_transcripts%>%rows_update(down)

  # get tables to create bar plots with

  Upregulated=table(all_transcripts$up)
  Downregulated=table(all_transcripts$down)


  full=rbind(Upregulated, Downregulated)

  # if Upregulated and Downregulated tables are not the same length, one must be extended to be the
  # same length as the other otherwise graphing doesn't work

  if(length(Upregulated)!=length(Downregulated)){
    if(length(Upregulated)<length(Downregulated)){
      diff=length(Downregulated)-length(Upregulated)
      for(i in 1:diff){
        full[[length(Upregulated)*2+1+(i-1)*2]]=0
      }
    }
  }

  # get largest bin from both upregulated & downregulated tables

  max=max(full)

  bp2=barplot(full, col=c("lightblue", "pink"), beside=TRUE,
              main="Upregulated & Downregulated Pause Sites by Transcript", xlab="Differential Pause Sites by Expression Change", ylab="Number of Transcripts", legend=c("Upregulated", "Downregulated"),
              args.legend=list(title="Pause Sites", pt.cex=1.3,cex=0.7,"topright", text.font=3),ylim=c(0,ceiling(max/100)*100))
  y=full[1:length(full)]

  # add labels for number of pause sites in each bin

  text(bp2,y+7,labels=as.character(y), cex=0.8)
}

#' Create bar plots to visualize the distribution of pause sites per transcript
#'
#' Takes as input a list of ribosome pause sites (the output of `get_pause_sites()`) and plots
#' a bar plot showing the distribution of pause sites, both constant and differential, per transcript.
#'
#' @param sites List of ribosome pause sites; usually the output of `get_pause_sites()`
#' @return `NULL`
#' @export
plot_ps_per_transcript=function(sites){
  sites=removezeros(sites)
  length=length(sites)

  ps=mcmapply(function(i){
    return(sites[[i]][[4]][[1]])
  }, i=1:length(sites), mc.cores=3)

  par(mfrow=c(1, 1))
  bp1=barplot(table(ps), col=c("gray"), main="Number of Differential Pause Sites per Transcript", xlab="Number of Differential Pause Sites", ylab="Number of Transcripts")
  text(bp1,0, table(ps),cex=0.9,pos=3)
}

#implement amino acid analysis

