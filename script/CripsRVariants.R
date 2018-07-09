library(CrispRVariants)
library(rtracklayer)
library(GenomicFeatures)
library(gridExtra)
txdb <- loadDb('./txdb/tair10_txdb.sqlite')

plates <- c('R1-R6_T2-1-1', 
           'R1-R6_T2-1-2',
           'R7-12_T1-1-1', 
           'R7-12_T1-1-2')
           
gd_fnames <- c('guide/gRNA_R1-R6.bed',
              'guide/gRNA_R1-R6.bed',
              'guide/gRNA_R19-R24.bed',
              'guide/gRNA_R19-R24.bed')

indel_table <- NULL

for (N in 1:length(plates)) {
  plate <- plates[N]
  gd_fname <- gd_fnames[N]
  
  gd <- rtracklayer::import(gd_fname)
  gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = 'center')
  
  for (i in 1:length(gdl)) {
    gname <- mcols(gdl)$name[i]
    folder <- paste0('figures/', plate, '/', gname)
    system(paste0('mkdir -p ', folder))
    for (RP in LETTERS[1:8]) {
      cat(paste0(plate, '_', gname, '_', RP, '...\n'))
      
      bam_fnames <- dir(paste0('bam/', plate), pattern = paste0(RP, '_[0-9]+.bam$'), full.names = TRUE)
      sample_names <- sub('.bam', '', dir(paste0('bam/', plate), pattern = paste0(RP, '_[0-9]+.bam$')))
      bam_fnames <- bam_fnames[c(1,5:12,2:4)]
      sample_names <- sample_names[c(1,5:12,2:4)]
      reference <- system(sprintf('samtools faidx index/tair10.fa %s:%s-%s', seqnames(gdl)[i], start(gdl)[i], end(gdl)[i]), intern = TRUE)[[2]]
      crispr_set <- readsToTarget(bam_fnames, target = gdl[i], reference = reference, names = sample_names, target.loc = 22, chimeras='ignore', verbose=FALSE)
      if (!is.null(crispr_set) && sum(grepl('no variant', rownames(crispr_set$cigar_freqs))) == 1) {
        pdf(paste0(folder, '/', plate, '_', gname, '_', RP, '.pdf'), wid=15, hei=15)
        p <- plotVariants(crispr_set, txdb=txdb, gene.text.size=8,
                          row.ht.ratio=c(1,8), col.wdth.ratio=c(4,2),
                          plotAlignments.args = list(line.weight=0.5, ins.size=2, legend.symbol.size=4),
                          plotFreqHeatmap.args = list(plot.text.size=3, x.size=8, legend.text.size=8,
                                                      legend.key.height=grid::unit(0.5, "lines")))
        dev.off()

        freqs <- crispr_set$cigar_freqs
        freqs_all <- colSums(freqs)
        freqs_indel <- colSums(freqs[-c(1, grep('SNV', rownames(freqs))),])
        indel_percent <- round(freqs_indel / freqs_all * 100, 2)
        sample_names <- names(freqs_all)
        indel_table <- rbind(indel_table, as.data.frame(cbind(plate, gname, sample_names, freqs_all, freqs_indel, indel_percent)))
      }
    }
  }
}

write.table(indel_table, './tables/indel_freq_all.tsv', row.names = F, quote=F, sep='\t')
