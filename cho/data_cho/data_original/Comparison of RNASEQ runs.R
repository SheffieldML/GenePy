#############################stadard starting lines####################
rm (list=ls())
####################################

#set file location 
setwd("E:/Dropbox/Shared folder-JLongworth & JGonzales/Data streams combined")



new_RNASEQ=read.delim("Galaxy27-[Cufflinks_on_data_20,_data_19,_and_data_25__transcript_expression].tabular")
old_RNASEQ=read.delim("Galaxy116-[Cufflinks_on_data_38,_data_10,_and_data_114__transcript_expression] (3).tabular")

new_RNASEQ2=new_RNASEQ[(new_RNASEQ$FPKM>0),]
old_RNASEQ2=old_RNASEQ[(old_RNASEQ$FPKM>0),]
data=merge.data.frame(new_RNASEQ2,old_RNASEQ2,by.y="tracking_id",by.x="tracking_id")

plot(data$FPKM.x,data$FPKM.y)


smoothScatter(data$FPKM.x,data$FPKM.y,ylim=c(0,6000),xlim=c(0,6000))

smoothScatter(data$FPKM.x,data$FPKM.y,ylim=c(0,60),xlim=c(0,60))
smoothScatter(data$FPKM.x,data$FPKM.y,ylim=c(0,10),xlim=c(0,10))
smoothScatter(log2(data$FPKM.x),log2(data$FPKM.y),ylim=c(-10,10),xlim=c(-10,10))


plot(data$FPKM.x,data$FPKM.y,ylim=c(0,1),xlim=c(0,1),pch=".")
cor(data$FPKM.x,data$FPKM.y)
