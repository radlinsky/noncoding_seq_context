all_miRNA<-"/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_GRCh37_primary_transcripts.BED"
table <- read.table(all_miRNA, header = F,stringsAsFactors = F,sep="\t")

# Very simply:
#  1) sort table by chr, position
#  2) for each chr:
#  3) go row by row, throw away an miRNA if it overlaps with the previous
chrs <- unique(table$V1)
keep_table <- data.frame()
toss_table <- data.frame()
for (chr in chrs){
  all_chr_rows <- which(table$V1 == chr)
  table_at_chr <- table[all_chr_rows,]
  table_at_chr <- table_at_chr[order(table_at_chr$V2),]
  start <- 0
  end <- 0
  keep <- c()
  toss <- c()
  for (row_i in seq(1,dim(table_at_chr)[1])){
    row <- table_at_chr[row_i,]
    if (row$V2 <= start || row$V2 <= end){
      toss <- c(toss, row_i)
      next
    }
    else {
      start = row$V2
    }
    if (row$V3 <= start || row$V3 <= end){
      toss <- c(toss, row_i)
      next
    }
    else {
      end = row$V3
      keep <- c(keep, row_i)
    }
    
  }
  keep_table <- rbind(keep_table, table_at_chr[keep,])
  toss_table <- rbind(toss_table, table_at_chr[toss,])
}

grouped_dir <- "/project/voight_subrate/cradens/noncoding_seq_context/data/results/mirbase/miRNA_grouped"
non_overlap_miRNA <- paste(grouped_dir,"/miRNA_GRCh37_primary_transcripts_non_overlap.BED",sep="")
write.table(keep_table, non_overlap_miRNA, row.names = F, col.names = F,quote = F,sep=",")