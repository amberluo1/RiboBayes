mmalias=mmrn$mmalias
mmsequences=data.frame()
for(i in 1:length(mmalias)){
  print(i)
  mmsequences=rbind(mmsequences, df%>%filter(grepl(mmalias[i], transcript)))
}
for(i in 1:nrow(eif5a)){
  mod=eif5a[[i,26]]
  if(mod==0){
    eif5a[[i,33]]=substr(eif5a[[i, 19]], floor(eif5a[[i,2]]/3)*3-30, floor(eif5a[[i,2]]/3)*3-28)
  }else if (mod==1){
    eif5a[[i,33]]=substr(eif5a[[i, 19]], floor(eif5a[[i,2]]/3)*3-32, floor(eif5a[[i,2]]/3)*3-30)
  }else{
    eif5a[[i,33]]=substr(eif5a[[i, 19]], floor(eif5a[[i,2]]/3)*3-31, floor(eif5a[[i,2]]/3)*3-29)
  }
}
for(i in 1:nrow(eif5a)){
  eif5a[[i,31]]=as.character(translate(DNAString(eif5a[[i,33]])))
                             }

eif5a%>%ggplot(mapping=aes(y=statistic))+geom_boxplot()+facet_wrap(~psite)
aachange=eif5a%>%group_by(psite)%>%summarize(stat=mean(statistic))
#create datasets for heatmap and proportion and zscore :D
eif5ad=eif5a%>%filter(statistic>0)
eif5au=eif5a%>%filter(statistic<0)

temp=as.data.frame(table(eif5a$rcodon))
temp=temp%>%rename(codon=Var1)%>%mutate(codon=as.character(codon))%>%rename(rnum=Freq)
codonusage2=codonusage2%>%left_join(temp, by="codon")
codonusage2[is.na(codonusage2)]=0
codonusage2=codonusage2%>%mutate(rsite=rnum/sum(rnum))

# codonusage=codonusage%>%left_join(temp, by="codon")
# codonusage[is.na(codonusage)]=0
# codonusage=codonusage%>%mutate(asite=anum/sum(anum))
aa2s=c()
pprop=aausage$psite
aprop=aausage$asite
bruh=c()
for(i in 1:21){
  for(j in 1:21){
    aa2=paste(aa[i], aa[j], sep="")
    proportion=pprop[i]*aprop[j]
    aa2s=c(aa2s, aa2)
    bruh=c(bruh, proportion)
  }
}
aa2usage=data.frame("aa2"=aa2s, "standard"=bruh)
aa2usage=aa2usage%>%left_join(codon2usage)
aa2usage[is.na(aa2usage)]=0

ptrans=c()
trans=eif5a$asite
for(i in 1:nrow(eif5a)){
  aa=as.character(translate(DNAString(eif5a[[i,29]])))
  ptrans=c(ptrans, aa)
}
ohno=which(trans!=ptrans)
for(i in ohno){
  mod=eif5a[[i,26]]
  if(mod==0){
    eif5a[[i,29]]=substr(eif5a[[i, 19]], floor(eif5a[[i,2]]/3)*3+15, floor(eif5a[[i,2]]/3)*3+17)
  }else if (mod==1){
    eif5a[[i,29]]=substr(eif5a[[i, 19]], floor(eif5a[[i,2]]/3)*3+13, floor(eif5a[[i,2]]/3)*3+15)
  }else{
    eif5a[[i,29]]=substr(eif5a[[i, 19]], floor(eif5a[[i,2]]/3)*3+14, floor(eif5a[[i,2]]/3)*3+16)
  }
}
glycine=eif5a%>%filter(acodon%in%c("GGT", "GGC", "GGA", "GGG"))
#clean for heatmap

for(i in 15){
  for(j in 1:nrow(codonusage2)){
    what=codonusage2[,i]
    what=what[which(is.finite(what))]
    min=min(what)
    what=codonusage2[,i]
    what=what[which(is.finite(what))]
    min=min(what)
    if(is.infinite(codonusage2[[j,i]])){
      codonusage2[[j,i]]=min
    }
  }
}
bruh=c("logp", "loge", "loga", "logr")
for(i in 7:10){
  for(j in 1:nrow(aausage2)){
    what=aausage2[[`i`]]
    what=what[which(is.finite(what))]
    min=min(what)
    if(is.infinite(aausage2[[j,i]])){
      aausage2[[j,i]]=min
    }
  }
}
#making heatmap

codonheatmap=codonusage2%>%select(logp, loge, loga, logr)
codons=codonusage2$codon
rownames(codonheatmap)=codons
codonheatmap=as.matrix(codonheatmap)
mydf=as.data.frame(codonusage$aa)
mydf=mydf%>%rename(`Amino Acid`=`codonusage$aa`)
rownames(mydf)=codonusage$codon
pheatmap(chm, col=redgreen(75), cexRow=0.5,fontsize_row=9,
         cexCol=1.25, cluster_cols = F, cluster_rows = F, breaks = seq(-5, 5, length.out = 76),gaps_row = c(2, 4, 6, 8, 14, 17, 21, 23, 25, 29, 32, 38, 39, 41, 45, 47, 53, 57, 61,62, 64
                                                                                              ),
         cellheight = 8, cellwidth = 25)

aausage2=as.data.frame(aausage2)
aminoheatmap=aausage2%>%select(logp, loge, loga, logr)
rownames(aminoheatmap)=aausage2$aa
aminoheatmap=as.matrix(aminoheatmap)
pheatmap(ahm2, color=redgreen(75), breaks = seq(-4.5, 4.5, length.out = 76),cexRow=0.5,fontsize_row=9,
         cexCol=1.25, cluster_cols = F, cluster_rows = F,
         cellheight = 12, cellwidth = 20)


heatmap.2(codonheatmap, Rowv=FALSE, Colv = FALSE, col=bluered(75), cexRow=0.6,breaks = seq(-4, 4, length.out = 100),
          cexCol=1.25, trace="none", density.info="none", key=FALSE, margins=c(4,3))

dumbbell(v1 = aa2usage$expected, v2 = aa2usage$actual, 
         text = TRUE, labels = aa2usage$aa2, segments = TRUE,
         pch = 19, pt.cex = 1.5, colv1 = "red", colv2="green")

logpdensity=data.frame("position"=c(-30:30), "negative","positive", "nonpolar")
for(j in 1:nrow(eif5a)){
  ack=rep(0, 61)
  for(i in -30:30){
  ack[i+31]=substr(eif5a[[j,21]], eif5a[[j, 22]]+i, eif5a[[j,22]]+i)
  }
  ack=as.data.frame(ack)
}
#density plots
yikes=data.frame("position"=c(-50:50), "positive"=0, "negative"=0, "ile"=0, "pro"=0, "gly"=0)
agh=data.frame()
for(i in 1:nrow(eif5a)){
  mini=data.frame("position"=c(-50:50), "aa"=0, "identifier"=eif5a[[i, 1]])
  for(j in -50:50){
    mini[[j+51, 2]]=substr(eif5a[[i, 21]], eif5a[[i, 22]]+j, eif5a[[i, 22]]+j)
  }
  agh=rbind(agh, mini)
}
positive=positive%>%mutate(position=as.numeric(as.character(Var1)), count=as.numeric(Freq), group="positive")%>%select(position, count, group)
ser=agh%>%filter(aa %in% c("S"))
ser=data.frame(table(ser$position))
ser=ser%>%mutate(position=as.numeric(as.character(Var1)), count=as.numeric(Freq), group="ser")%>%select(position, count, group)
ser=ser%>%mutate(counts=count/731, standard=0.080, logFC=log(counts/standard, base=2), label=NA)
density=rbind(density, ser)
ggplot(density%>%filter(abs(position)<25), aes(x = position, y = logFC, color=group)) +geom_area(alpha=0.8,aes(fill=group))+
      geom_line()+facet_wrap(~group)+theme_minimal()+
      theme(legend.position="none")+geom_vline(xintercept=0, linetype="dashed")
