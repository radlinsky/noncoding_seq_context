# The file path of the miRNA file:
mirna_file = "/project/voight_subrate/cradens/noncoding_seq_context/data/not_mine/hsa.gff3"
out_file = "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh37.BED"

mirna<-read.table(mirna_file,sep="\t",skip=13,header=F,stringsAsFactors=F)

# remove y-chromosome miRNA
mirna = subset(mirna,mirna$V1!="chrY")

# remove column 2 if it is empty
if (length(unique(mirna$V2)) == 1){
  mirna$V2<-NULL
}

# remove column 6 if it is empty
if (length(unique(mirna$V6)) == 1){
  mirna$V6<-NULL
}

# remove column 8 if it is empty
if (length(unique(mirna$V8)) == 1){
  mirna$V8<-NULL
}

# re-order columns:
# chr | start | end | typr | strand | name/id
mirna <- mirna[,c(1,3,4,2,5,6)]

write.table(mirna,out_file,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")