cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggplot2)
library(forcats)
library(zoo)

give.n <- function(x){
  return(c(y = 1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


# True data EBOLA
reftrue = read.table("edges_ref.txt",header=T,sep="\t",na.strings="N/A")
reftrueweight = read.table("edges_ref_weight.txt",header=T,sep="\t",na.strings="N/A")
reftruenocollapse = read.table("edges_ref_nocollapse.txt",header=T,sep="\t",na.strings="N/A")
refphyrax = read.table("comptreesrax.txt",header=T,sep="\t",na.strings="N/A")
refphyraxnocoll = read.table("comptreesrax_nocollapse.txt",header=T,sep="\t",na.strings="N/A")

refphyrax=refphyrax[refphyrax$found=="true" & refphyrax$topodepth>1,]
refphyraxnocoll=refphyraxnocoll[refphyraxnocoll$found=="true" & refphyraxnocoll$topodepth>1,]
refphyrax$collapse="true"
refphyraxnocoll$collapse="false"

cor(refphyrax$support,refphyrax$comparedsupport)
cor(refphyrax$support,refphyrax$comparedsupport,method="spearman")
cor(refphyraxnocoll$support,refphyraxnocoll$comparedsupport)
cor(refphyraxnocoll$support,refphyraxnocoll$comparedsupport,method="spearman")

refphyrax=rbind(refphyrax,refphyraxnocoll)

svg("Correlation_phyml_raxml_ebola.svg",width=8,height=4.5)
ggplot(refphyrax,aes(x=support,y=comparedsupport))+geom_point(size=0.7)+theme_bw()+facet_wrap(~collapse)
dev.off()


len=18996

reftrue[order(reftrue$length),"sortidx"] = seq(1,length(reftrue$length)) 
reftrueweight[order(reftrueweight$length),"sortidx"] = seq(1,length(reftrueweight$length)) 
reftruenocollapse[order(reftruenocollapse$length),"sortidx"] = seq(1,length(reftruenocollapse$length)) 

svg("ebola_panel_brlen_internal.svg",width=4.7,height=2.2)
tmplen=reftrue[reftrue$terminal=="false", ]
tmplen[order(tmplen$length),"sortidx"]=seq(1,length(tmplen$length))/length(tmplen$length)
ggplot(tmplen,aes(y=sortidx,x=len*length))+geom_point(size=0.10)+theme_bw()+xlim(c(0,10))+scale_x_continuous(limits=c(0,10),breaks=0:10)
dev.off()

theo_values=data.frame(nmutstr=c("1","2","1","2"),support=c(0.632,0.865,0.905,0.995),boot=c("Frequentist","Frequentist","Bayesian","Bayesian"))
theo_values$nmutstr=factor(theo_values$nmutstr,levels = c("0","1","2",">=2"))

reftrueweight$boot="Bayesian"
reftrue$boot="Frequentist"
reftruenocollapse$boot="FrequentistNoColl"
reftrue2=rbind(reftrue,reftrueweight,reftruenocollapse)
reftrue2$nmuts=round(reftrue2$length*len)
reftrue2$nmutstr[reftrue2$nmuts==0]="0"
reftrue2$nmutstr[reftrue2$nmuts==1]="1"
reftrue2$nmutstr[reftrue2$nmuts==2]="2"
reftrue2$nmutstr[reftrue2$nmuts>2]=">=2"
reftrue2$nmutstr=factor(reftrue2$nmutstr,levels = c("0","1","2",">=2"))
reftrue2=reftrue2[reftrue2$boot=="Bayesian" | reftrue2$boot=="Frequentist" | (reftrue2$boot=="FrequentistNoColl" & reftrue2$nmuts==0),]
reftrue2=reftrue2[reftrue2$topodepth>1,]
reftrue2sup0=reftrue2[reftrue2$nmuts>0,]

svg("ebola_panel_1.svg",width=8,height=3)
ggplot(reftrue2, aes(x=factor(boot),y=support,fill=boot)) + 
  geom_boxplot(outlier.size = 0.5) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  #geom_point(data=reftrue2sup0,aes(factor(boot),y=support),alpha=0.4,color="black",pch = 21, position = position_jitterdodge(jitter.width=0.2,dodge.width = 0.75))+
  facet_grid(cols=vars(nmutstr)) + 
  geom_hline(data=theo_values, aes(yintercept=support,color=boot))+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()
dev.off()


c1=reftrue2[reftrue2$terminal=="false" & reftrue2$boot=="Frequentist","support"]
c2=reftrue2[reftrue2$terminal=="false" & reftrue2$boot=="Bayesian","support"]
print(cor(c1,c2))
print(cor(c1,c2,method = "spearman"))

