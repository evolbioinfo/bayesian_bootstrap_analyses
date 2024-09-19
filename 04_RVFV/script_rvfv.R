cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggplot2)
library(forcats)
library(zoo)

give.n <- function(x){
  return(c(y = 1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# True data RVFV
reftrueS = read.table("results_S/edges_ref.txt",header=T,sep="\t",na.strings="N/A")
reftrueweightS = read.table("results_S/edges_ref_weight.txt",header=T,sep="\t",na.strings="N/A")
reftruenocollapseS = read.table("results_S/edges_ref_nocollapse.txt",header=T,sep="\t",na.strings="N/A")
reftrueS$alilen=1736 # S
reftrueweightS$alilen=1736 # S
reftruenocollapseS$alilen=1736 # S

reftrueM = read.table("results_M/edges_ref.txt",header=T,sep="\t",na.strings="N/A")
reftrueweightM = read.table("results_M/edges_ref_weight.txt",header=T,sep="\t",na.strings="N/A")
reftruenocollapseM = read.table("results_M/edges_ref_nocollapse.txt",header=T,sep="\t",na.strings="N/A")
reftrueM$alilen=3967
reftrueweightM$alilen=3967
reftruenocollapseM$alilen=3967

reftrueL = read.table("results_L/edges_ref.txt",header=T,sep="\t",na.strings="N/A")
reftrueweightL = read.table("results_L/edges_ref_weight.txt",header=T,sep="\t",na.strings="N/A")
reftruenocollapseL = read.table("results_L/edges_ref_nocollapse.txt",header=T,sep="\t",na.strings="N/A")
reftrueL$alilen=6471
reftrueweightL$alilen=6471
reftruenocollapseL$alilen=6471

reftrue=rbind(reftrueS, reftrueM,reftrueL)
reftrueweight=rbind(reftrueweightS,reftrueweightM,reftrueweightL)
reftruenocollapse=rbind(reftruenocollapseS,reftruenocollapseM,reftruenocollapseL)

reftrue[order(reftrue$length),"sortidx"] = seq(1,length(reftrue$length)) 
reftrueweight[order(reftrueweight$length),"sortidx"] = seq(1,length(reftrueweight$length)) 
reftruenocollapse[order(reftruenocollapse$length),"sortidx"] = seq(1,length(reftruenocollapse$length)) 

svg("rvfv_brlen.svg",width=4.7,height=2.2)
tmplen=reftrue[reftrue$terminal=="false", ]
tmplen[order(tmplen$length*tmplen$alilen),"sortidx"]=seq(1,length(tmplen$length))/length(tmplen$length)
ggplot(tmplen,aes(y=sortidx,x=alilen*length))+geom_point(size=0.10)+theme_bw()+xlim(c(0,10))+scale_x_continuous(limits=c(0,10),breaks=0:10)
dev.off()

theo_values=data.frame(nmutstr=c("1 mutation","2 mutations","1 mutation","2 mutations"),support=c(0.632,0.865,0.905,0.995),boot=c("Frequentist","Frequentist","Bayesian","Bayesian"))
theo_values$nmutstr=factor(theo_values$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))

reftrueweight$boot="Bayesian"
reftrue$boot="Frequentist"
reftruenocollapse$boot="NoCollapse"
reftrue2=rbind(reftrue,reftrueweight,reftruenocollapse)
reftrue2$nmuts=round(reftrue2$length*reftrue2$alilen)
reftrue2$nmutstr[reftrue2$nmuts==0]="0 mutation"
reftrue2$nmutstr[reftrue2$nmuts==1]="1 mutation"
reftrue2$nmutstr[reftrue2$nmuts==2]="2 mutations"
reftrue2$nmutstr[reftrue2$nmuts>2]=">=3 mutations"
reftrue2$nmutstr=factor(reftrue2$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))
reftrue2=reftrue2[reftrue2$boot=="Bayesian" | reftrue2$boot=="Frequentist" | reftrue2$boot=="NoCollapse",]
reftrue2=reftrue2[reftrue2$topodepth>1,]
reftrue2sup0=reftrue2[reftrue2$nmuts>0,]

svg("rvfv_supports.svg",width=7.7,height=3.2)
ggplot(reftrue2, aes(x=factor(boot),y=support,fill=boot)) + 
  geom_boxplot(outlier.size = 0.5) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  facet_grid(cols=vars(nmutstr)) + 
  geom_hline(data=theo_values, aes(yintercept=support,color=boot))+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()
dev.off()

c1=reftrue2[reftrue2$terminal=="false" & reftrue2$boot=="Frequentist" & reftrue2$nmuts>0,"support"]
c2=reftrue2[reftrue2$terminal=="false" & reftrue2$boot=="Bayesian" & reftrue2$nmuts>0,"support"]
c3=reftrue2[reftrue2$terminal=="false" & reftrue2$boot=="NoCollapse" & reftrue2$nmuts>0,"support"]
print(cor(c1,c2))
print(cor(c3,c2))


