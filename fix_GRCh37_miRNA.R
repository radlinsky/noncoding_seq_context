########
# LiftOver throws away some column data outside of chr|start|end for some reason.
# This script recovers those columns by going back to the original BED file and
# identifying which rows were not unMapped by liftOver, and gets the column data.

mirna37 <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh37.BED"
mirna37 <- read.table(mirna37,sep="\t",header=F,stringsAsFactors=F)
mirna38 <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh38.BED"
mirna38 <- read.table(mirna38,sep="\t",header=F,stringsAsFactors=F)
original <- "/project/voight_subrate/cradens/noncoding_seq_context/data/not_mine/hsa.gff3"
original <- read.table(original,sep="\t",skip=13,header=F,stringsAsFactors=F)

# remove y-chromosome miRNA
original = subset(original,original$V1!="chrY")

# remove column 2 if it is empty
if (length(unique(original$V2)) == 1){
  original$V2<-NULL
}

# remove column 6 if it is empty
if (length(unique(original$V6)) == 1){
  original$V6<-NULL
}

# remove column 8 if it is empty
if (length(unique(original$V8)) == 1){
  original$V8<-NULL
}

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

mirna37$type <- rep("EMPTY",length(mirna37$V1))
mirna37$strand <- rep("EMPTY",length(mirna37$V1))

for (id in not_lost_in_maplation){
  # what row in original file contained non-lost id?
  id_index <- which(original$V9 == id)
  # what was the original row's corresponding type and strand?
  orig_type <- original$V3[id_index]
  orig_strand <- original$V7[id_index]
  
  id_index <- which(mirna37$V4 == id)
  mirna37$type[id_index] <- orig_type
  
  
  mirna37$strand[id_index] <- orig_strand
}

if ("EMPTY" %in% mirna37$type){
  stop("'EMPTY' found in type column. This is not expected.")
}
if ("EMPTY" %in% mirna37$strand){
  stop("'EMPTY' found in strand column. This is not expected.")
}

fixed_mirna <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh37_full.BED"
write.table(x=mirna37,
            file=fixed_mirna,
            quote=FALSE,
            sep="\t",
            row.names=FALSE,
            col.names = FALSE)

# mature miRNA = "MIMA####"
mature_mirna<- grepl("MIMA",mirna37$V4)
# stem_loop precursors = "MI#####"
primary_transcripts <- mirna37[!mature_mirna,]

mirna_primary_transcripts <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh37_primary_transcripts.BED"
write.table(x=primary_transcripts,
            file=mirna_primary_transcripts,
            quote=FALSE,
            sep="\t",
            row.names=FALSE,
            col.names = FALSE)