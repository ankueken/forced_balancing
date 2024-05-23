args=(commandArgs(TRUE))
library(R.matlab)
data=readMat(args)

S=data$model[,,1][[1]]
rev=rep(0,dim(S)[2])
dim(S)

source('deficiency')
library(igraph)

D=deficiency(S,rev)

write.table(D$A,file=paste(substring(args,1,nchar(args)-4),'.A',sep=''))
write.table(D$Y,file=paste(substring(args,1,nchar(args)-4),'.Y',sep=''))
