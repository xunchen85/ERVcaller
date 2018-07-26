#main <- function() {
 args <- commandArgs(trailingOnly = TRUE)
 filename <- args[1]
  dat <- read.delim(file = filename, header = FALSE,sep="")
  insert_size <- dat
  insert_size <- abs(subset(insert_size,insert_size[,1]=="P")[,13])
  Quantile<-quantile(insert_size,0.1)
  Quantile<-data.frame(Quantile)
  dat[,2]<-Quantile$Quantile
 write.table(dat,file=filename,sep="\t",col.names = F, row.names = F,quote = FALSE);
#}
