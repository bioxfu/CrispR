files <- dir('bam2', pattern = '*.stats', recursive = TRUE, full.names = TRUE)
sampleID <- sub('.bam.stats', '', dir('bam2', pattern = '*.stats', recursive = TRUE))

counts <- NULL

for(i in 1:length(files)) {
  dfm <- read.table(files[i])[c(1,3)]
  dfm$Sample <- sampleID[i]
  counts <- rbind(counts, dfm[c(3,1,2)])
}

colnames(counts) <- c('sample', 'gRNA', 'count')
counts$gRNA <- sub('_.+', '', counts$gRNA)

write.table(counts, 'tables/gRNA_counts.tsv', sep = '\t', row.names = FALSE)
