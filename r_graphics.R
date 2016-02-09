
table <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/ensemble_grch37/miRNA/summary_SNVs.txt"
table<-read.table(table,stringsAsFactors = F, header = T, sep=",")
options(bitmapType='cairo')
for (pop_i in 2:length(names(table))){
  pop = names(table)[pop_i]
  plot_name = paste("/project/voight_subrate/cradens/noncoding_seq_context/data/results/ensemble_grch37/miRNA/",pop, "_miRNA_SNVs_Observed.png",sep="")
  png(plot_name, width=10, height=7.5, units="in", res=200)
  zero_data=table[,pop_i][which(table[,pop_i] ==0)]
  non_zero_data=table[,pop_i][which(table[,pop_i] !=0)]
  hist(non_zero_data,main=paste("miRNA SNVs for",pop,"(",length(zero_data),"0s)"),xlim=c(0,10),xlab="Number of SNVs per miRNA")
  dev.off()
}



