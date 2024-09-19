library(concom)

args = commandArgs(trailingOnly=TRUE)

infile=args[1]
outfile=args[2]

# Read distance matrix
dist=read.table(infile,skip=1,header=F,row.names = 1)
colnames(dist)=row.names(dist)

# Convert distances in adjacency matrix
# distances of 0 converted to 1 (link: identical)
# distances >0 converted to 0 (no link: not identical)
adjMat=dist
adjMat[adjMat==0]=-1
adjMat[adjMat>0]=0
adjMat[adjMat==-1]=1

# Compute connected components
comp=concomFromMatAdj(as.matrix(adjMat))

# We keep the first sequences of each component
tokeep=unlist(lapply(1:length(comp$components),function(c){unlist(comp$components[[c]][[1]])}))
namestokeep=row.names(adjMat)[tokeep]

# We write the list of sequence names to keep in the end
write.table(namestokeep,file=outfile, col.names = F, row.names = F, quote = FALSE)
