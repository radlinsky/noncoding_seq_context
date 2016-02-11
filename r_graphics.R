
table <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/ensemble_grch37/lincRNA/summary_SNVs.txt"
table<-read.table(table,stringsAsFactors = F, header = T, sep=",")
options(bitmapType='cairo')
for (pop_i in 2:length(names(table))){
  pop = names(table)[pop_i]
  plot_name = paste("/project/voight_subrate/cradens/noncoding_seq_context/data/results/ensemble_grch37/lincRNA/",pop, "_lincRNA_SNVs_Observed.png",sep="")
  png(plot_name, width=10, height=7.5, units="in", res=200)
  zero_data=table[,pop_i][which(table[,pop_i] ==0)]
  non_zero_data=table[,pop_i][which(table[,pop_i] !=0)]
  hist(table[,pop_i],main=paste("lincRNA SNVs for",pop),xlim=c(0,10000),xlab="Number of SNVs per lincRNA")
  dev.off()
}



png("plot_name.png", width=10, height=7.5, units="in", res=200)
nt->
dev.off()

for (pop_i in 2:length(names(table))){
  pop = names(table)[pop_i]
  plot_name = paste("/project/voight_subrate/cradens/noncoding_seq_context/data/results/ensemble_grch37/lincRNA/",pop, "_box_lincRNA_SNVs_Observed.png",sep="")
  png(plot_name, width=10, height=7.5, units="in", res=200)
  b
  dev.off()
}