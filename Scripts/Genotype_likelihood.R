#main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  filename <- args[1]
  filename2 <- args[2]
  array1 <- read.delim(file = filename, header = FALSE,sep="")
#  for (i in 1:nrow(array1)) {
    array1[,12]<-as.numeric(as.character(array1[,10]))
    array1[,13]<-as.numeric(as.character(array1[,10]))+as.numeric(as.character(array1[,11]))
    k <-array1[,12]
    n <-array1[,13]
    array1[,14]<-round((dbinom(array1[,12],array1[,13],0.95))/(dbinom(array1[,12],array1[,13],0.95)+dbinom(array1[,12],array1[,13],0.5)+dbinom(array1[,13]-array1[,12],array1[,13],0.95)),4)
    array1[,15]<-round((dbinom(array1[,12],array1[,13],0.5))/(dbinom(array1[,12],array1[,13],0.95)+dbinom(array1[,12],array1[,13],0.5)+dbinom(array1[,13]-array1[,12],array1[,13],0.95)),4)
    array1[,16]<-round((dbinom(array1[,13]-array1[,12],array1[,13],0.95))/(dbinom(array1[,12],array1[,13],0.95)+dbinom(array1[,12],array1[,13],0.5)+dbinom(array1[,13]-array1[,12],array1[,13],0.95)),4)
    array1[,17]<-round(10*(-log10(1-array1[,14])))
    array1[,17]<-ifelse (array1[,17]>=40,40,array1[,17])
    array1[,18]<-round(10*(-log10(1-array1[,15])))
    array1[,18]<-ifelse (array1[,18]>=40,40,array1[,18])
    array1[,19]<-round(10*(-log10(1-array1[,16])))
    array1[,19]<-ifelse (array1[,19]>=40,40,array1[,19])
#  }
  write.table(array1,file=filename2,sep="\t",row.names = F,col.names = F,quote = FALSE);
#}
