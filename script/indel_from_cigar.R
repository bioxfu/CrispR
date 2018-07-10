dfm <- read.table('tables/reads_from_gRNA_with_cigar', stringsAsFactors = F)
colnames(dfm) <- c('plate', 'gname', 'sample_names', 'cigar')
g <- paste(dfm$plate, dfm$gname, dfm$sample_names, sep = '\t')

indel_freq <- function(x) {
  freqs_all <- round(length(x), 0)
  freqs_indel <- round(length(grep('[ID]', x)), 0)
  indel_percent <- round(freqs_indel / freqs_all * 100, 2)
  c(freqs_all, freqs_indel, indel_percent)
}
result <- do.call(rbind, tapply(X = dfm$cigar, INDEX = g, indel_freq))
result <- data.frame(rownames(result), result)
colnames(result) <- c('plate\tgname\tsample_names', 'freqs_all', 'freqs_indel', 'indel_percent')

write.table(result, 'tables/indel_freq_all_double_check.tsv', row.names = F, quote = F, sep = '\t')
