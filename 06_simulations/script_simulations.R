cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggplot2)
library(zoo)

give.n <- function(x){
  return(c(y = 1, label = length(x))) 
}

give.true <- function(x){
  return(c(y = -0.1, label = length(x))) 
}


fullrefscalelen=data.frame()
fullrefweightscalelen=data.frame()
fullrefnocollapsescalelen=data.frame()
for(scale in c(1,4,16,64)){
  len=floor(18996/scale)

  reflen = read.table(paste0("results_",scale,"_",len,"/edges_ref_refcol.txt"),header=T,sep="\t",na.strings="N/A")
  refweightlen = read.table(paste0("results_",scale,"_",len,"/edges_ref_weight_refcol.txt"),header=T,sep="\t",na.strings="N/A")
  refnocollapselen = read.table(paste0("results_",scale,"_",len,"/edges_ref_nocollapse_refcol.txt"),header=T,sep="\t",na.strings="N/A")

  reflen$scale=scale
  refweightlen$scale=scale
  refweightlen[order(refweightlen$length),"sortidx"]=seq(1,length(refweightlen$length))
  reflen$len=len
  reflen[order(reflen$length),"sortidx"]=seq(1,length(reflen$length))
  refweightlen$len=len
  refnocollapselen$scale=scale
  refnocollapselen$len=len
  refnocollapselen[order(refnocollapselen$length),"sortidx"]=seq(1,length(refnocollapselen$length))
  
  fullrefscalelen=rbind(fullrefscalelen,reflen)
  fullrefweightscalelen=rbind(fullrefweightscalelen,refweightlen)
  fullrefnocollapsescalelen=rbind(fullrefnocollapsescalelen,refnocollapselen)
}

fullrefweightscalelen$boot="Bayesian"
fullrefscalelen$boot="Frequentist"
fullrefnocollapsescalelen$boot="NoCollapse"

fullrefscalelen2=rbind(fullrefscalelen,fullrefweightscalelen,fullrefnocollapsescalelen)
fullrefscalelen2$nmuts=round(fullrefscalelen2$length*fullrefscalelen2$len)
fullrefscalelen2$nmutstr[fullrefscalelen2$nmuts==0]="0 mutation"
fullrefscalelen2$nmutstr[fullrefscalelen2$nmuts==1]="1 mutation"
fullrefscalelen2$nmutstr[fullrefscalelen2$nmuts==2]="2 mutations"
fullrefscalelen2$nmutstr[fullrefscalelen2$nmuts>2]=">=3 mutations"
fullrefscalelen2$nmutstr=factor(fullrefscalelen2$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))
fullrefscalelen2=fullrefscalelen2[fullrefscalelen2$boot=="Bayesian" | fullrefscalelen2$boot=="Frequentist" | (fullrefscalelen2$boot=="NoCollapse" & fullrefscalelen2$nmuts==0),]
fullrefscalelen2=fullrefscalelen2[fullrefscalelen2$topodepth>1,]
fullrefscalelen2$found[fullrefscalelen2$nmuts==0]="true"
fullrefscalelen2found=fullrefscalelen2[fullrefscalelen2$found=="false" & fullrefscalelen2$nmuts>0,]

theo_values=data.frame(nmutstr=c("1 mutation","2 mutations","1 mutation","2 mutations"),support=c(0.632,0.865,0.905,0.995),boot=c("Frequentist","Frequentist","Bayesian","Bayesian"))
theo_values$nmutstr=factor(theo_values$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))

svg("simulated_panel_1.svg",width=10,height=7)
ggplot(fullrefscalelen2, aes(x=factor(scale),y=support,fill=boot, alpha=found)) +
  geom_boxplot(outlier.size = 0.1,color="black") + 
  #geom_jitter(width = 0.2)+
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  stat_summary(data=fullrefscalelen2[fullrefscalelen2$found=="true",], fun.data = give.true, geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  facet_wrap(~nmutstr) + 
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  scale_alpha_discrete(range = c(0.3, 1.0))+
  geom_hline(data=theo_values, aes(yintercept=support,fill=boot,color=boot))+
  scale_y_continuous(limits=c(-0.1,1),breaks=seq(0,1,0.1))+
  theme_bw()
dev.off()

# Correlations per simulations
for(f in c(1,4,16,64)){
  c1=fullrefscalelen2[fullrefscalelen2$terminal=="false" & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Frequentist","support"]
  c2=fullrefscalelen2[fullrefscalelen2$terminal=="false" & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Bayesian","support"]
  print(paste(f,cor(c1,c2)))
  print(paste(f,cor(c1,c2,method = "spearman")))
}

# AUC per simulation
for(f in c(1,4,16,64)){
  # FREQ
  onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Frequentist",c("found","support")]
  TP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support>=x,"support"]))}))
  TPTOTAL=length(onemut[onemut$found=="true","support"])
  FP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support>=x,"support"]))}))
  FPTOTAL=length(onemut[onemut$found=="false","support"])
  plot(FP,TP, type='l',col="blue")
  print(paste0("Freq, F=",f," AUC=",sprintf("%.3f",sum(diff(rev(FP/FPTOTAL))*rollmean(rev(TP/TPTOTAL),2)))))
  
  #BAYES
  onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Bayesian",c("found","support")]
  TP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support>=x,"support"]))}))
  TPTOTAL=length(onemut[onemut$found=="true","support"])
  FP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support>=x,"support"]))}))
  FPTOTAL=length(onemut[onemut$found=="false","support"])
  lines(FP,TP)
  print(paste0("Bayes, F=",f," AUC=",sprintf("%.3f",sum(diff(rev(FP/FPTOTAL))*rollmean(rev(TP/TPTOTAL),2)))))
  #readline()
}

# AUC for union of all simulations
# FREQ
fullrefscalelen3=rbind(fullrefscalelen,fullrefweightscalelen,fullrefnocollapsescalelen)
fullrefscalelen3$nmuts=round(fullrefscalelen3$length*fullrefscalelen3$len)
fullrefscalelen3$nmutstr[fullrefscalelen3$nmuts==0]="0 mutation"
fullrefscalelen3$nmutstr[fullrefscalelen3$nmuts==1]="1 mutation"
fullrefscalelen3$nmutstr[fullrefscalelen3$nmuts==2]="2 mutations"
fullrefscalelen3$nmutstr[fullrefscalelen3$nmuts>2]=">=3 mutations"
fullrefscalelen3$nmutstr=factor(fullrefscalelen3$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))
fullrefscalelen3=fullrefscalelen3[fullrefscalelen3$boot=="Bayesian" | fullrefscalelen3$boot=="Frequentist" | fullrefscalelen3$boot=="NoCollapse",]
fullrefscalelen3=fullrefscalelen3[fullrefscalelen3$topodepth>1,]
fullrefscalelen3$found[fullrefscalelen3$nmuts==0]="true"
fullrefscalelen3found=fullrefscalelen3[fullrefscalelen3$found=="false" & fullrefscalelen3$nmuts>0,]

svg("Simu_roc.svg",width=5,height=5.5)
mutcutoff=1
onemut=fullrefscalelen3[fullrefscalelen3$nmuts>=mutcutoff & fullrefscalelen3$boot=="Frequentist",c("found","support")]
TP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support>=x,"support"]))}))
TPTOTAL=length(onemut[onemut$found=="true","support"])
FP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support>=x,"support"]))}))
TN=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support<x,"support"]))}))
FN=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support<x,"support"]))}))
FPTOTAL=length(onemut[onemut$found=="false","support"])
type1=FP/(TN+FP)
power=TP/(TP+FN)
data.frame(cut=seq(0,1,0.01),tp=TP,fp=FP,type1=type1,power=power)
plot(FP,TP, type='l',col="blue",lwd = 2)
print(paste0("Freq, F=ALL"," AUC=",sprintf("%.3f",sum(diff(rev(FP/FPTOTAL))*rollmean(rev(TP/TPTOTAL),2)))))
#BAYES
onemut=fullrefscalelen3[fullrefscalelen3$nmuts>=mutcutoff & fullrefscalelen3$boot=="Bayesian",c("found","support")]
TP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support>=x,"support"]))}))
TPTOTAL=length(onemut[onemut$found=="true","support"])
FP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support>=x,"support"]))}))
TN=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support<x,"support"]))}))
FN=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support<x,"support"]))}))
FPTOTAL=length(onemut[onemut$found=="false","support"])
type1=FP/(TN+FP)
power=TP/(TP+FN)
data.frame(cut=seq(0,1,0.01),tp=TP,fp=FP,type1=type1,power=power)
lines(FP,TP,lwd = 2)
print(paste0("Bayes, F=ALL"," AUC=",sprintf("%.3f",sum(diff(rev(FP/FPTOTAL))*rollmean(rev(TP/TPTOTAL),2)))))
#NOCOLLAPSE
onemut=fullrefscalelen3[fullrefscalelen3$nmuts>=mutcutoff & fullrefscalelen3$boot=="NoCollapse",c("found","support")]
TP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support>=x,"support"]))}))
TPTOTAL=length(onemut[onemut$found=="true","support"])
FP=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support>=x,"support"]))}))
TN=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="false" & onemut$support<x,"support"]))}))
FN=unlist(lapply(seq(0,1,0.01),function(x){return(length(onemut[onemut$found=="true" & onemut$support<x,"support"]))}))
FPTOTAL=length(onemut[onemut$found=="false","support"])
type1=FP/(TN+FP)
power=TP/(TP+FN)
data.frame(cut=seq(0,1,0.01),tp=TP,fp=FP,type1=type1,power=power)
lines(FP,TP,lwd = 2, col="orange")
print(paste0("NoCollapse, F=ALL"," AUC=",sprintf("%.3f",sum(diff(rev(FP/FPTOTAL))*rollmean(rev(TP/TPTOTAL),2)))))
dev.off()


# Type1 & Recall
for(f in c(1,4,16,64)){
  # Cutoff
  for(c in c(0.5,0.631,0.7,0.9)){
    # FREQ
    onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Frequentist",c("found","support")]
    TP=length(onemut[onemut$found=="true" & onemut$support>=c,"support"])
    FP=length(onemut[onemut$found=="false" & onemut$support>=c,"support"])
    FN=length(onemut[onemut$found=="true" & onemut$support<c,"support"])
    TN=length(onemut[onemut$found=="false" & onemut$support<c,"support"])
    type1=FP/(TN+FP)
    power=TP/(TP+FN)
    print(paste0("Freq:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Type1:",sprintf("%.3f",type1)))
    print(paste0("Freq:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Recall:",sprintf("%.3f",power)))
    #BAYES
    onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Bayesian",c("found","support")]
    # Cutoff
    TP=length(onemut[onemut$found=="true" & onemut$support>=c,"support"])
    FP=length(onemut[onemut$found=="false" & onemut$support>=c,"support"])
    FN=length(onemut[onemut$found=="true" & onemut$support<c,"support"])
    TN=length(onemut[onemut$found=="false" & onemut$support<c,"support"])
    type1=FP/(TN+FP)
    power=TP/(TP+FN)
    print(paste0("Bayes:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Type1:",sprintf("%.3f",type1)))
    print(paste0("Bayes:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Recall:",sprintf("%.3f",power)))
  }
  #readline()
}

# Screen for Type1~0.1 per simulation
for(f in c(1,4,16,64)){
  # Cutoff
  for(c in seq(0,1,0.01)){
    # FREQ
    onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Frequentist",c("found","support")]
    TP=length(onemut[onemut$found=="true" & onemut$support>=c,"support"])
    FP=length(onemut[onemut$found=="false" & onemut$support>=c,"support"])
    FN=length(onemut[onemut$found=="true" & onemut$support<c,"support"])
    TN=length(onemut[onemut$found=="false" & onemut$support<c,"support"])
    type1=FP/(TN+FP)
    power=TP/(TP+FN)
    if(type1>=0.04 && type1<=0.06){
      print(paste0("Freq:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Type1:",sprintf("%.3f",type1)))      
      print(paste0("Freq:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Recall:",sprintf("%.3f",power)))
    }

    #BAYES
    onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$scale==f & fullrefscalelen2$boot=="Bayesian",c("found","support")]
    # Cutoff
    TP=length(onemut[onemut$found=="true" & onemut$support>=c,"support"])
    FP=length(onemut[onemut$found=="false" & onemut$support>=c,"support"])
    FN=length(onemut[onemut$found=="true" & onemut$support<c,"support"])
    TN=length(onemut[onemut$found=="false" & onemut$support<c,"support"])
    type1=FP/(TN+FP)
    power=TP/(TP+FN)
    if(type1>=0.04 && type1<=0.06){
      print(paste0("Bayes:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Type1:",sprintf("%.3f",type1)))      
      print(paste0("Bayes:",f," cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Recall:",sprintf("%.3f",power)))
    }
  }
}

# Screen for Type1~0.1 for union of all simulations
for(c in seq(0,1,0.01)){
  # FREQ
  onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$boot=="Frequentist",c("found","support")]
  TP=length(onemut[onemut$found=="true" & onemut$support>=c,"support"])
  FP=length(onemut[onemut$found=="false" & onemut$support>=c,"support"])
  FN=length(onemut[onemut$found=="true" & onemut$support<c,"support"])
  TN=length(onemut[onemut$found=="false" & onemut$support<c,"support"])
  type1=FP/(TN+FP)
  power=TP/(TP+FN)
  if(type1>=0.04 && type1<=0.06){
    print(paste0("Freq: cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Type1:",sprintf("%.3f",type1)))      
    print(paste0("Freq: cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Recall:",sprintf("%.3f",power)))
  }
  
  #BAYES
  onemut=fullrefscalelen2[fullrefscalelen2$nmuts>=1 & fullrefscalelen2$boot=="Bayesian",c("found","support")]
  # Cutoff
  TP=length(onemut[onemut$found=="true" & onemut$support>=c,"support"])
  FP=length(onemut[onemut$found=="false" & onemut$support>=c,"support"])
  FN=length(onemut[onemut$found=="true" & onemut$support<c,"support"])
  TN=length(onemut[onemut$found=="false" & onemut$support<c,"support"])
  type1=FP/(TN+FP)
  power=TP/(TP+FN)
  if(type1>=0.04 && type1<=0.06){
    print(paste0("Bayes: cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Type1:",sprintf("%.3f",type1)))      
    print(paste0("Bayes: cut:",c," ","TP:",TP,"|FP:",FP,"|TN:",TN,"|FN:",FN," Recall:",sprintf("%.3f",power)))
  }
}


## Branch lengths
svg("simulated_panel_2_internal.svg",width=5,height=4)
tmplen1=fullrefscalelen[fullrefscalelen$scale==1 & fullrefscalelen$terminal=="false", ]
tmplen64=fullrefscalelen[fullrefscalelen$scale==64 & fullrefscalelen$terminal=="false", ]
tmplen1[order(tmplen1$length),"sortidx"]=seq(1,length(tmplen1$length))/length(tmplen1$length)
tmplen64[order(tmplen64$length),"sortidx"]=seq(1,length(tmplen64$length))/length(tmplen64$length)
tmplen=rbind(tmplen1,tmplen64)
ggplot(tmplen,aes(y=sortidx,x=len*length))+geom_point(size=0.10)+facet_grid(scale~.)+theme_bw()+xlim(c(0,10))+scale_x_continuous(limits=c(0,10),breaks=0:10)
dev.off()

