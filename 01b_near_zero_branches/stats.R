library(ggplot2)
library(ggthemes)

args = commandArgs(trailingOnly=TRUE)

indir=args[1]

setwd(paste0(indir))

# Import supports statistics of phyml collapsed & non collapsed
phymlnocollapsesupports=read.table("phyml_nocollapse_supports.txt",header=T,sep="\t",na.strings = "N/A")
phymlcollapsesupports=read.table("phyml_collapse_supports.txt",header=T,sep="\t",na.strings = "N/A")

# Import supports statistics of iq-tree collapsed & non collapsed
iqtreenocollapsesupports=read.table("iqtree_nocollapse_supports.txt",header=T,sep="\t",na.strings = "N/A")
iqtreecollapsesupports=read.table("iqtree_collapse_supports.txt",header=T,sep="\t",na.strings = "N/A")

# Import supports statistics of raxml collapsed & non collapsed
raxmlnocollapsesupports=read.table("raxml_nocollapse_supports.txt",header=T,sep="\t",na.strings = "N/A")
raxmlcollapsesupports=read.table("raxml_collapse_supports.txt",header=T,sep="\t",na.strings = "N/A")

# Import branch comparison between iq-tree & phyml, collapsed & non collapsed
iqphynocollapsesupports=read.table("phyml_iqtree_nocollapse_supports.txt",header=T,sep="\t",na.strings = "N/A")
iqphynocollapsesupports=iqphynocollapsesupports[iqphynocollapsesupports$found=="true" & iqphynocollapsesupports$topodepth>1,]
iqphycollapsesupports=read.table("phyml_iqtree_collapse_supports.txt",header=T,sep="\t",na.strings = "N/A")
iqphycollapsesupports=iqphycollapsesupports[iqphycollapsesupports$found=="true" & iqphycollapsesupports$topodepth>1,]

# Import branch comparison between raxml & phyml, collapsed & non collapsed
raxphynocollapsesupports=read.table("phyml_raxml_nocollapse_supports.txt",header=T,sep="\t",na.strings = "N/A")
raxphynocollapsesupports=raxphynocollapsesupports[raxphynocollapsesupports$found=="true" & raxphynocollapsesupports$topodepth>1,]
raxphycollapsesupports=read.table("phyml_raxml_collapse_supports.txt",header=T,sep="\t",na.strings = "N/A")
raxphycollapsesupports=raxphycollapsesupports[raxphycollapsesupports$found=="true" & raxphycollapsesupports$topodepth>1,]


# We bind all the rows, for support stats
tmpsupports=phymlnocollapsesupports[,c("length","support")]
tmpsupports$tool="PhyML"
tmpsupports$collapse=FALSE
rawsupports=tmpsupports

tmpsupports=phymlcollapsesupports[,c("length","support")]
tmpsupports$tool="PhyML"
tmpsupports$collapse=TRUE
rawsupports=rbind(rawsupports,tmpsupports)

rawsupports$support=rawsupports$support*100 # Same scale as IQTREE

tmpsupports=iqtreenocollapsesupports[,c("length","support")]
tmpsupports$tool="IQTREE"
tmpsupports$collapse=FALSE
rawsupports=rbind(rawsupports,tmpsupports)

tmpsupports=iqtreecollapsesupports[,c("length","support")]
tmpsupports$tool="IQTREE"
tmpsupports$collapse=TRUE
rawsupports=rbind(rawsupports,tmpsupports)


tmpsupports=raxmlnocollapsesupports[,c("length","support")]
tmpsupports$tool="RAXML"
tmpsupports$collapse=FALSE
tmpsupports$support=tmpsupports$support*100
rawsupports=rbind(rawsupports,tmpsupports)

tmpsupports=raxmlcollapsesupports[,c("length","support")]
tmpsupports$tool="RAXML"
tmpsupports$collapse=TRUE
tmpsupports$support=tmpsupports$support*100
rawsupports=rbind(rawsupports,tmpsupports)

ggplot(rawsupports,aes(x=length*18994,y=support))+geom_point(size=0.7)+facet_wrap(tool ~ collapse,ncol = 2)+theme_bw()+xlim(0,10)+scale_x_continuous(limits = c(0, 10), breaks = 0:10)

rawsupports$NullBranch="Non_Null_Branches"
rawsupports$NullBranch[rawsupports$length<=5.264820469621986e-06]="Null_Branches"

rawsupports$NullBranch=factor(rawsupports$NullBranch,levels = c("Null_Branches","Non_Null_Branches"))

svg("null_branches.svg",width=7,height=3)
ggplot(rawsupports,aes(y=support,x=tool,fill=tool,alpha=collapse))+geom_violin(draw_quantiles = c(0.5))+facet_wrap(~NullBranch)+scale_colour_colorblind()+theme_bw()
dev.off()

table(rawsupports[!is.na(rawsupports$support),c("collapse","tool","NullBranch")])
table(rawsupports[!is.na(rawsupports$support) & rawsupports$support>=70,c("collapse","tool","NullBranch")])


# We bind the rows, for IQTREE vs PhyML  branch comparisons
tmpsupport=iqphynocollapsesupports[,c("support","comparedsupport")]
tmpsupport$collapse=FALSE
colnames(tmpsupport)=c("PhyML_Support","IQTREE_Support","Collapse")
compsupports = tmpsupport
tmpsupport=iqphycollapsesupports[,c("support","comparedsupport")]
tmpsupport$collapse=TRUE
colnames(tmpsupport)=c("PhyML_Support","IQTREE_Support","Collapse")
compsupports = rbind(compsupports,tmpsupport)

compsupports$PhyML_Support=compsupports$PhyML_Support*100

ggplot(compsupports,aes(x=PhyML_Support,y=IQTREE_Support))+geom_point(size=0.7)+facet_wrap(~ Collapse)+theme_bw()

cor(iqphynocollapsesupports$support,iqphynocollapsesupports$comparedsupport,use="complete.obs")
cor(iqphycollapsesupports$support,iqphycollapsesupports$comparedsupport,use="complete.obs")


# We bind the rows, for RAxML vs PhyML  branch comparisons
tmpsupport=raxphynocollapsesupports[,c("support","comparedsupport")]
tmpsupport$collapse=FALSE
colnames(tmpsupport)=c("PhyML_Support","RAxML_Support","Collapse")
compsupports = tmpsupport
tmpsupport=raxphycollapsesupports[,c("support","comparedsupport")]
tmpsupport$collapse=TRUE
colnames(tmpsupport)=c("PhyML_Support","RAxML_Support","Collapse")
compsupports = rbind(compsupports,tmpsupport)


compsupports$PhyML_Support=compsupports$PhyML_Support*100
compsupports$RAxML_Support=compsupports$RAxML_Support*100

ggplot(compsupports,aes(x=PhyML_Support,y=RAxML_Support))+geom_point(size=0.7)+facet_wrap(~ Collapse)+theme_bw()

cor(raxphynocollapsesupports$support,raxphynocollapsesupports$comparedsupport,use="complete.obs")
cor(raxphycollapsesupports$support,raxphycollapsesupports$comparedsupport,use="complete.obs")
