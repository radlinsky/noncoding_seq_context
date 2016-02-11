########
# LiftOver throws away some column data outside of chr|start|end for some reason.
# This script recovers those columns by going back to the original BED file and
# identifying which rows were not unMapped by liftOver, and gets the column data.

mirna37 <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh37.BED"
mirna37 <- read.table(mirna37,sep="\t",header=F,stringsAsFactors=F)
mirna38 <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh38.BED"
mirna38 <- read.table(mirna38,sep="\t",header=F,stringsAsFactors=F)
original < - "/project/voight_subrate/cradens/noncoding_seq_context/data/not_mine/hsa.gff3"
original <- read.table(original,sep="\t",skip=13,header=F,stringsAsFactors=F)

not_lost_in_maplation <- c()
lost_in_maplation <- c()
for(id in mirna38[,4]){
  if (id %in% mirna37[,4]){
    not_lost_in_maplation <- c(not_lost_in_maplation, id)
  }
  # Else the row was thrown away by liftOver for whatever reason
  else{
    lost_in_maplation <- c(lost_in_maplation, id)
  }
}


