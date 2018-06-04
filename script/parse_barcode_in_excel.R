library(readxl)
library(magrittr)

input_files <- sub('.xlsx', '', dir('tables', '*.xlsx', full.names = T))

for (file in input_files) {
  
  dfm <- read_excel(paste0(file, '.xlsx'), sheet = 3, skip = 1)
  
  FP_primer <- dfm[grep('FP-', dfm[, 1, drop=T]),]
  RP_primer <- dfm[grep('RP-', dfm[, 1, drop=T]),]
  colnames(FP_primer) <- c('primer', 'FP_sequence')
  colnames(RP_primer) <- c('primer', 'RP_sequence')
  
  FP_primer$primer <- sub('FP-', '', FP_primer$primer)
  RP_primer$primer <- sub('RP-', '', RP_primer$primer)
  
  FP_primer$FP_sequence <- sub("5'-ACTCTTTCCCTACACGACGCTCTTCCGATCTgctt", '', FP_primer$FP_sequence)
  FP_primer$FP_sequence <- sub("tggagtgagtacggtgtgc", '', FP_primer$FP_sequence)
  RP_primer$RP_sequence <- sub("5'-GACTGGAGTTCAGACGTGTGCTCTTCCGATCTctgt", '', RP_primer$RP_sequence)
  RP_primer$RP_sequence <- sub("tgagttggatgctggatgg", '', RP_primer$RP_sequence)
  
  sample_all <- data.frame(RP=rep(RP_primer$primer, each=12), 
                           FP=rep(FP_primer$primer, 8), stringsAsFactors = F)
  sample_all <- merge(sample_all, FP_primer, by.x = 2, by.y = 1)
  sample_all <- merge(sample_all, RP_primer, by.x = 2, by.y = 1)
  sample_all$FP <- as.numeric(sample_all$FP)
  sample_all <- sample_all[order(sample_all$RP, sample_all$FP), c('RP', 'FP', 'RP_sequence', 'FP_sequence')]
  
  write.table(sample_all, paste0(file, '.barcode_sequence.tsv'), sep='\t', quote=F, row.names = F)

}
