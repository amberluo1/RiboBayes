zscorestats=mcmapply(function(x,y,z,w){
  return(zscorestats_lp(x,y,z,w))
}, x = pkmtranscripts[1:30], y = size, z=threshold, w=ssize,
mc.cores=7)

zscorepeakstats=function(transcript, size, threshold, ssize){
  cov=getpcov(transcript)
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=laggedtransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  set=makeset(acov, acov2, threshold)
  positions=set$position
  count=length(positions)
  pcov=cov%>%filter(position %in% positions)
  set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
                 "control1"=1:count, "control2"=1:count, "control3"=1:count)
  set[,1]=positions
  for(i in 1:count){
    for(j in 2:7){
      set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)]
    }
  }
  set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
    mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  set2 <- set[,-1]
  rownames(set2) <- set[,1]
  kappaset=set
  ah=differentialpeaks(set)
  diffpeaks=ah[[2]]
  diffpositions=diffpeaks[,c(1:3)]
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
  allpeaks=c(rep(0,9))
  diffpeaks=c(rep(0,9))
  kappaset=kappaset%>%select(-avgKD, -avgcontrol, -ratio)
  kappaset=kappaset%>%mutate(kd=0, control=0)
  
  if(!(nrow(kappaset)==0)){
    for(i in 1:nrow(kappaset)){
      for(j in 1:2){
        for(k in 1:3){
          if(kappaset[[i, (j-1)*3+k+1]]==1){
            kappaset[[i, 7+j]]=kappaset[[i, 7+j]]+1
          }
        }
      }
    }
    for(i in 1:nrow(kappaset)){
      if((kappaset[[i, 8]]+kappaset[[i, 9]])==1){
        allpeaks[[1]]=allpeaks[[1]]+1
      }
      else if((kappaset[[i, 8]]==1 & kappaset[[i, 9]]==1)){
        allpeaks[[2]]=allpeaks[[2]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==0)){
        allpeaks[[3]]=allpeaks[[3]]+1
      }
      else if((kappaset[[i, 8]]==1 & kappaset[[i, 9]]==2)){
        allpeaks[[4]]=allpeaks[[4]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==0)){
        allpeaks[[5]]=allpeaks[[5]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==2)){
        allpeaks[[6]]=allpeaks[[6]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==1)){
        allpeaks[[7]]=allpeaks[[7]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==2)){
        allpeaks[[8]]=allpeaks[[8]]+1
      }
      else if((kappaset[[i, 8]]==0 & kappaset[[i, 9]]==2)){
        allpeaks[[9]]=allpeaks[[9]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==1)){
        allpeaks[[10]]=allpeaks[[10]]+1
      }
      else if((kappaset[[i, 8]]==0 & kappaset[[i, 9]]==3)){
        allpeaks[[11]]=allpeaks[[11]]+1
      }
      else if((kappaset[[i, 8]]==1 & kappaset[[i, 9]]==3)){
        allpeaks[[12]]=allpeaks[[12]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==3)){
        allpeaks[[13]]=allpeaks[[13]]+1
      }
      else{
        allpeaks[[14]]=allpeaks[[14]]+1
      }
    }
  }
  
  ahh=diffpositions
  diffpositions=diffpositions$position
  diffkappaset=kappaset%>%filter(position %in% diffpositions)
  diffkappaset=diffkappaset%>%mutate(kd=0, control=0)
  
  if(!(length(diffpositions)==0)){
    for(i in 1:nrow(diffkappaset)){
      for(j in 1:2){
        for(k in 1:3){
          if(diffkappaset[[i, (j-1)*3+k+1]]==1){
            diffkappaset[[i, 7+j]]=diffkappaset[[i, 7+j]]+1
          }
        }
      }
    }
    for(i in 1:nrow(diffkappaset)){
      if((diffkappaset[[i, 8]]+diffkappaset[[i, 9]])==1){
        diffpeaks[[1]]=diffpeaks[[1]]+1
      }
      else if((diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==1)){
        diffpeaks[[2]]=diffpeaks[[2]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==0)){
        diffpeaks[[3]]=diffpeaks[[3]]+1
      }
      else if((diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[4]]=diffpeaks[[4]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==0)){
        diffpeaks[[5]]=diffpeaks[[5]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[6]]=diffpeaks[[6]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==1)){
        diffpeaks[[7]]=diffpeaks[[7]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[8]]=diffpeaks[[8]]+1
      }
      else if((diffkappaset[[i, 8]]==0 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[9]]=diffpeaks[[9]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==1)){
        diffpeaks[[10]]=diffpeaks[[10]]+1
      }
      else if((diffkappaset[[i, 8]]==0 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[11]]=diffpeaks[[11]]+1
      }
      else if((diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[12]]=diffpeaks[[12]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[13]]=diffpeaks[[13]]+1
      }
      else{
        diffpeaks[[14]]=diffpeaks[[14]]+1
      }
    }
  }
  allpos=kappaset$position
  numpeaks=c(nrow(set), nrow(diffkappaset))
  newlist=list(allpeaks, diffpeaks, numpeaks, allpos, ahh)
  return(newlist)
}

zscorestats_op=function(transcript, size, threshold, ssize){
  cov=getpcov(transcript)
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=zscoretransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  return(get_icc(set))
}
zscorestats_lp=function(transcript, size, threshold, ssize){
  cov=getpcov(transcript)
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=laggedtransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  if(nrow(set)<=2){
    return(0)
  }
  return(get_icc(set))
}

getzscorestats_l=function(transcript, size, threshold, ssize){
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 29,
                      range.upper = 36,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = get_experiments(pkm.ribo))
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=laggedtransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  if(nrow(set)<=2){
    return(0)
  }
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  return(get_icc(set))
}

getzscorestats_o=function(transcript, size, threshold, ssize){
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 29,
                      range.upper = 36,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = get_experiments(pkm.ribo))
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=zscoretransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  return(get_icc(set))
}

laggedtransform=function(cov, window, threshold, subwindow){
  acov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    arrange(experiment)%>%mutate(scan1=0, scan2=0, scan3=0)
  max=acov%>%arrange(desc(position))
  max=max[[1,2]]
  ohboy=acov%>%arrange(desc(position))
  for(i in 1:6){
    print(i)
    values=acov$counts
    max=ohboy[[1,2]]
    max2=ohboy[[1,2]]
    j3=window
    for(j in window:max){
      position=(i-1)*max2+j3
      position2=(i-1)*max2+j
      mean=mean(values[(position-window+1):(position-4)])
      sd=sd(values[(position-window+1):(position-4)])
      acov[[position2, 5]]=(values[position]-mean)/sd
      if(is.na(acov[[position2, 5]] | is.infinite(acov[[position2, 5]]))){
        acov[[position2, 5]]=0
      }
      if(acov[[position2, 5]]>threshold & values[position]>7){
        print(j)
        values=values[-position]
        max=max-1
        j3=j3-1
      }
      j3=j3+1
    }
  }
  
  for(i in 1:6){
    values=acov$counts
    max=ohboy[[1,2]]
    max2=ohboy[[1,2]]
    print(i)
    for(j in (max-window+1):1){
      position=(i-1)*max2+j
      mean=mean(values[(position+4):(position+3+window)])
      sd=sd(values[(position+4):(position+3+window)])
      acov[[position, 6]]=(values[position]-mean)/sd
      if(is.na(acov[[position, 6]] | is.infinite(acov[[position, 6]]))){
        acov[[position, 6]]=0
      }
      if(acov[[position, 6]]>threshold & values[position]>threshold){
        print(j)
        values=values[-position]
      }
    }
  }
  for(i in 1:6){
    values=acov$counts
    max=ohboy[[1,2]]
    print(i)
    for(j in (max-window/2+1):(window/2)){
      position=(i-1)*max+j
      mean=mean(values[c((position-window/2+1):(position-3), (position+3):(position+window/2-1))])
      sd=sd(values[c((position-window/2+1):(position-3), (position+3):(position+window/2-1))])
      acov[[position, 7]]=(values[position]-mean)/sd
      if(is.na(acov[[position, 7]]) | is.infinite(acov[[position, 7]])){
        acov[[position, 7]]=0
      }
    }
  }
  acov=acov%>%mutate(z=(scan1+scan2+scan3)/3)
  acov2=acov%>%filter(z>threshold & counts>7)
  acov=acov[,c(1:6, 8, 7)]
  return(list(acov, acov2))
}
zscoretransform=function(cov, size, threshold, ssize){
  acov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    mutate(avgcount=mean(count))%>%
    arrange(experiment)%>%mutate(countr=counts+0.0001, z=counts)
  acov2=acov%>%arrange(desc(position))
  max=acov2[[1,2]]
  #calculates z - leaves out the beginning and end of each transcript
  window=c(1:size)
  subwindow=c(1:ssize)
  stuff=array(dim=c(4, max, 6))
  matrix=array(dim=c(3, 3))
  matrix[,1]=c(1, 0, 0)
  matrix[,2]=c(0, 1, -1)
  matrix[,3]=c(0, 1, 0)
  
  for(i in 0:5){ #experiment number
    print(i)
    for(j in 0:2){ #moving window position
      f=matrix[1, j+1]
      g=matrix[2, j+1]
      h=matrix[3, j+1]
      for(k in (j*(size/2)+f):(max-size+j*(size/2)+f)){#moving window range of nucleotide positions
        m=1
        save=c(1:2*ssize+1)
        counting=1
        for(l in (-1*j*size/2+g):((size-1-j*(size/2))+g+h)){#mini-loop to calculate mean and SD
          window[m]=as.numeric(acov[[i*max+k+l, 6]])
          m=m+1
        }
        if(j==0){
          window=window[-c(1:ssize)]
        }
        else if(j==1){
          window=window[-c((size/2-ssize+2):(size/2+ssize))]
        }
        else{
          window=window[-c((size-ssize+1):(size))]
        }
        mean=mean(window)
        sd=sd(window)
        zscore=(acov[[i*max+k, 7]]-mean)/(sd+0.0001)
        stuff[j+1,k,i+1]=zscore
      }
    }
  }
  acov2=acov%>%arrange(desc(position))
  max=acov2[[i,2]]
  values=c(1:3)
  
  #calculating 4th row of 3D array
  for(j in 1:6){
    for(i in 1+size:max-size){
      for(k in 1:3){
        values[k]=as.numeric(stuff[k, i, j])
      }
      mean=mean(values)
      if(is.na(mean)){
        mean=0
      }
    }
  }
  #setting last column of acov
  for(i in 0:5){
    for(j in (size):(max-size)){
      stuff[4, j, i+1]=(stuff[1, j, i+1]+stuff[2, j, i+1]+stuff[3, j, i+1])/3
      acov[i*max+j, 7]=stuff[4, j, i+1]
    }
  }
  #acov=acov%>%mutate(score=z*counts)
  acov2=acov%>%filter(z>6 & counts>7)
  newlist = list(acov, acov2)
  return (newlist)
}


positions=set$position
count=length(positions)
cov_peaks=cov%>%filter(position %in% positions)
set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
               "control1"=1:count, "control2"=1:count, "control3"=1:count)
set[,1]=positions
for(i in 1:count){
  for(j in 2:7){
    set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)-1]
  }
}
set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
  mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))%>%mutate(sig=0)
set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
m=1
for(i in 1:nrow(set)){
  if(set[[m,8]]>5 | set[[m,9]]>5){
    set[[m, 11]]=1
  }
  if(set[[m,11]]==0){
    set=set[-m,]
    m=m-1
  }
  m=m+1
}
set=set%>%select(-sig)
set2 <- set[,-1]
rownames(set2) <- set[,1]
positions=set$position
count=length(positions)
pcov=cov%>%filter(position %in% positions)
set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
               "control1"=1:count, "control2"=1:count, "control3"=1:count)
set[,1]=positions
for(i in 1:count){
  for(j in 2:7){
    set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)]
  }
}
set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
  mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))
set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
