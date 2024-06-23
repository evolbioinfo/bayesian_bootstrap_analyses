cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggplot2)
library(forcats)
library(zoo)

give.n <- function(x){
  return(c(y = 1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

give.true <- function(x){
  return(c(y = -0.1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

fullrefscalelen=data.frame()
fullrefweightscalelen=data.frame()
fullrefnocollapsescalelen=data.frame()
for(scale in c(1,4,16,64)){
  len=floor(18996/scale)
  reflen = read.table(paste0("results_",scale,"_",len,"/edges_ref.txt"),header=T,sep="\t",na.strings="N/A")
  refweightlen = read.table(paste0("results_",scale,"_",len,"/edges_ref_weight.txt"),header=T,sep="\t",na.strings="N/A")
  refnocollapselen = read.table(paste0("results_",scale,"_",len,"/edges_ref_nocollapse.txt"),header=T,sep="\t",na.strings="N/A")

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
fullrefscalelen2=fullrefscalelen2[fullrefscalelen2$topodepth>1,]
fullrefscalelen2=fullrefscalelen2[fullrefscalelen2$nmuts==0,]

theo_values=data.frame(nmutstr=c("1 mutation","2 mutations","1 mutation","2 mutations"),support=c(0.632,0.865,0.905,0.995),boot=c("Frequentist","Frequentist","Bayesian","Bayesian"))
theo_values$nmutstr=factor(theo_values$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))

svg("Supp_Fig_2.svg",width=7,height=7)
ggplot(fullrefscalelen2, aes(x=boot,y=support,fill=boot, alpha=found)) +
  geom_boxplot(outlier.size = 0.1,color="black") + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  stat_summary(data=fullrefscalelen2[fullrefscalelen2$found=="true",], fun.data = give.true, geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  facet_wrap(~scale) + 
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  scale_alpha_discrete(range = c(0.3, 1.0))+
  geom_hline(data=theo_values, aes(yintercept=support,fill=boot,color=boot))+
  scale_y_continuous(limits=c(-0.1,1),breaks=seq(0,1,0.1))+
  theme_bw()
dev.off()
