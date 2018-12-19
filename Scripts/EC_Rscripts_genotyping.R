#main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  filename <- args[1]
  filename2 <- args[2]
  dat <- read.delim(file = filename, header = FALSE,sep="")
  dat <- subset(dat,dat[,2]>0)
  list <- read.delim(file = filename2, header = TRUE,sep="")
  Quantile<-quantile(dat[,2]/2,0.2)
  Quantile<-data.frame(Quantile)
  list[,4]<-Quantile$Quantile
  list$SD<-sd(dat[,2]/2)
  list$Num<-nrow(dat)
  list$Mean<-mean(dat[,2]/2)
  write.table(list,file=filename2,sep="\t",row.names = F,quote = FALSE);
#}
