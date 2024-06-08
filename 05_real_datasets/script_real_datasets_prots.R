cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

library(ggplot2)
library(forcats)
library(zoo)
library(dplyr)
library(plyr)

give.n <- function(x){
  return(c(y = 1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

samples=read.table("homoplasies_prot.txt",header=F)
colnames(samples)=c("name","alignpars","phymlpars","parshomo","mlhomo")

fullref=data.frame()
fullrefweight=data.frame()
fullrefnocollapses=data.frame()
for(s in samples$name){
  print(s)
  reftrue = read.table(paste0(s,"/edges_ref.txt"),header=T,sep="\t",na.strings="N/A")
  reftrueweight = read.table(paste0(s,"/edges_ref_weight.txt"),header=T,sep="\t",na.strings="N/A")
  reftruenocollapse = read.table(paste0(s,"/edges_ref_nocollapse.txt"),header=T,sep="\t",na.strings="N/A")
  alilen=read.table(paste0(s,"/alilen.txt"),header=F,sep="\t",na.strings="N/A")
  ntax=read.table(paste0(s,"/ntaxa.txt"),header=F,sep="\t",na.strings="N/A")
  alphabet=read.table(paste0(s,"/alphabet.txt"),header=F,sep="\t",na.strings="N/A")
  homoplasy=read.table(paste0(s,"/homoplasy.txt"),header=F,sep="\t",na.strings="N/A")
  mlhomoplasy=read.table(paste0(s,"/mlhomoplasy.txt"),header=F,sep="\t",na.strings="N/A")
  
  ntax=ntax[1,1]
  alilen=alilen[1,1]
  alphabet=alphabet[1,1]
  homoplasy=homoplasy[1,1]
  mlhomoplasy=mlhomoplasy[1,1]
  
  print(ntax)
  print(alilen)
  print(alphabet)
  
  reftrue$sample=s
  reftrueweight$sample=s
  reftruenocollapse$sample=s
  
  reftrue$alphabet=alphabet
  reftrueweight$alphabet=alphabet
  reftruenocollapse$alphabet=alphabet
  
  reftrue$homoplasy=homoplasy
  reftrueweight$homoplasy=homoplasy
  reftruenocollapse$homoplasy=homoplasy

  reftrue$mlhomoplasy=mlhomoplasy
  reftrueweight$mlhomoplasy=mlhomoplasy
  reftruenocollapse$mlhomoplasy=mlhomoplasy
    
  reftrue$alilen=alilen
  reftrueweight$alilen=alilen
  reftruenocollapse$alilen=alilen

  reftrue$ntaxa=ntax
  reftrueweight$ntaxa=ntax
  reftruenocollapse$ntaxa=ntax
  
  fullref=rbind(fullref,reftrue)
  fullrefweight=rbind(fullrefweight,reftrueweight)
  fullrefnocollapses=rbind(fullrefnocollapses,reftruenocollapse)
}

tmplenhomo=fullref[fullref$terminal=="false", ]
tmplenhomo[order(tmplenhomo$length*tmplenhomo$alilen),"sortidx"]=seq(1,length(tmplenhomo$length))
tmplenhomo$sortidx=tmplenhomo$sortidx/length(tmplenhomo$sortidx)
svg("Real_datasets_brlen_protein_all.svg",width=4.7,height=2.2)
ggplot(tmplenhomo,aes(y=sortidx,x=alilen*length))+geom_point(size=0.10)+theme_bw()+xlim(c(0,10))+scale_x_continuous(limits=c(0,10),breaks=0:10)
dev.off()

theo_values=data.frame(nmutstr=c("1 mutation","2 mutations","1 mutation","2 mutations"),support=c(0.632,0.865,0.905,0.995),boot=c("Frequentist","Frequentist","Bayesian","Bayesian"))
theo_values$nmutstr=factor(theo_values$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))

fullrefweight$boot="Bayesian"
fullref$boot="Frequentist"
fullrefnocollapses$boot="FrequentistNoColl"

# All datasets
reftrue2=rbind(fullref,fullrefweight,fullrefnocollapses)
reftrue2$nmuts=round(reftrue2$length*reftrue2$alilen)
reftrue2$nmutstr[reftrue2$nmuts==0]="0 mutation"
reftrue2$nmutstr[reftrue2$nmuts==1]="1 mutation"
reftrue2$nmutstr[reftrue2$nmuts==2]="2 mutations"
reftrue2$nmutstr[reftrue2$nmuts>2]=">=3 mutations"
reftrue2$nmutstr=factor(reftrue2$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))
reftrue2=reftrue2[reftrue2$boot=="Bayesian" | reftrue2$boot=="Frequentist" | (reftrue2$boot=="FrequentistNoColl" & reftrue2$nmuts==0),]
reftrue2=reftrue2[reftrue2$topodepth>1,]
reftrue2sup0=reftrue2[reftrue2$nmuts>0,]
svg("Real_datasets_supports_all.svg",width=7.7,height=3)
ggplot(reftrue2, aes(x=factor(boot),y=support,fill=boot)) + 
  geom_boxplot(outlier.size = 0.5) + 
  stat_summary(fun.data = give.n, geom = "text", fun = median, position = position_dodge(width = 0.75)) +
  facet_grid(cols=vars(nmutstr)) + 
  geom_hline(data=theo_values, aes(yintercept=support,color=boot))+
  scale_fill_manual(values=cbbPalette)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()
dev.off()

tmpref = fullref
tmpref$supportweight= fullrefweight$support
tmpref=tmpref[tmpref$topodepth>1,]
tmpref$nmuts=round(tmpref$length*tmpref$alilen)
tmpref$nmutstr[tmpref$nmuts==0]="0 mutation"
tmpref$nmutstr[tmpref$nmuts==1]="1 mutation"
tmpref$nmutstr[tmpref$nmuts==2]="2 mutations"
tmpref$nmutstr[tmpref$nmuts>2]=">=3 mutations"
tmpref$nmutstr=factor(tmpref$nmutstr,levels = c("0 mutation","1 mutation","2 mutations",">=3 mutations"))
cors <- ddply(tmpref, c("nmutstr","alphabet"), summarise, cor = round(cor(support, supportweight), 2))
corspear <- ddply(tmpref, c("nmutstr","alphabet"), summarise, cor = round(cor(support, supportweight,method="spearman"), 2))

c1=tmpref[tmpref$terminal=="false","support"]
c2=tmpref[tmpref$terminal=="false","supportweight"]
print(cor(c1,c2))
print(cor(c1,c2,method = "spearman"))
