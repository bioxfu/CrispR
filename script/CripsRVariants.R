library(CrispRVariants)
library(rtracklayer)
library(GenomicFeatures)
library(gridExtra)

argv <- commandArgs(T)
txdb <- loadDb(argv[1])
fasta <- argv[2]
gRNAtable <- read.table(argv[3], stringsAsFactors = F)
output_indel <- argv[4]
output_snv <- argv[5]
output_DI <- argv[6]
output_aln <- argv[7]
pam <- argv[8]

# txdb <- loadDb('~/Gmatic7/gene/tair10/txdb/tair10_txdb.sqlite')
# fasta <- '~/Gmatic7/genome/tair10/tair10.fa'
# gRNAtable <- read.table('guide/gRNA.tsv', stringsAsFactors = F)
# output_indel <- 'test_indel.tsv'
# output_snv <- 'test_snv'
# output_DI <- 'test_DI'
# output_aln <- 'test_aln'
# pam <- 'Cas9'

plates <- gRNAtable$V1
gd_fnames <- gRNAtable$V2

indel_table <- NULL

for (N in 1:length(plates)) {
  plate <- plates[N]
  gd_fname <- gd_fnames[N]
  
  gd <- rtracklayer::import(gd_fname)
  gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = 'center')
  
  tot <- data.frame(matrix(nrow=96, ncol=2))
  rownames(tot) <- paste0(rep(LETTERS[1:8], each=12), '_', 1:12)
  tot2 <- read.table(paste0('split/', plate, '/reads_stat.tsv'), row.names = 1)
  tot <- merge(tot, tot2, by.x=0, by.y=0, all.x=T)
  rownames(tot) <- tot[,1]
  tot <- tot[paste0(rep(LETTERS[1:8], each=12), '_', 1:12), c(4,5)]
  colnames(tot) <- c('barcode', 'total_reads')
  
  for (i in 1:length(gdl)) {
    gname <- mcols(gdl)$name[i]
    folder <- paste0('figures/', plate, '/', gname)
    system(paste0('mkdir -p ', folder))

    snv_table <- as.data.frame(matrix(0, nrow = 4*29, ncol = 96))
    colnames(snv_table) <- paste0(rep(LETTERS[1:8], each=12), '_', 1:12)
    rownames(snv_table) <- paste0(rep(c(-23:-1, 1:6), each=4), c('A', 'T', 'C', 'G'))
    DI_table <- data.frame(NULL)
    
    alns <- rep(0, 96)
    names(alns) <- paste0(rep(LETTERS[1:8], each=12), '_', 1:12)
    for (RP in LETTERS[1:8]) {
      cat(paste0(plate, '_', gname, '_', RP, '...\n'))
      
      bam_fnames <- dir(paste0('bam/', plate), pattern = paste0(RP, '_[0-9]+.bam$'), full.names = TRUE)
      sample_names <- sub('.bam', '', dir(paste0('bam/', plate), pattern = paste0(RP, '_[0-9]+.bam$')))
      bam_fnames <- bam_fnames[c(1,5:12,2:4)]
      sample_names <- sample_names[c(1,5:12,2:4)]
      #reference <- system(sprintf('/home/xfu/miniconda2/envs/gmatic/bin/samtools faidx %s %s:%s-%s', fasta, seqnames(gdl)[i], start(gdl)[i], end(gdl)[i]), intern = TRUE)[[2]]
      reference <- system(sprintf('samtools faidx %s %s:%s-%s', fasta, seqnames(gdl)[i], start(gdl)[i], end(gdl)[i]), intern = TRUE)[[2]]
      if (as.character(strand(gdl[i])) == '-') {
        reference <- Biostrings::reverseComplement(Biostrings::DNAString(reference))
      }
      if (pam == 'Cas9') {
        # target.loc = 5 + 17
        crispr_set <- readsToTarget(bam_fnames, target = gdl[i], reference = reference, names = sample_names, target.loc = 22, chimeras='ignore', verbose=FALSE)
      }
      if (pam == 'Cpf1') {
        # target.loc = 5 + (4 + 19)
        crispr_set <- readsToTarget(bam_fnames, target = gdl[i], reference = reference, names = sample_names, target.loc = 28, chimeras='ignore', verbose=FALSE)
      }
      
      if (!is.null(crispr_set)) {
        pdf(paste0(folder, '/', plate, '_', gname, '_', RP, '.pdf'), wid=15, hei=15)
        if (pam == 'Cas9') {
          p <- plotVariants(crispr_set, txdb=txdb, gene.text.size=8,
                            row.ht.ratio=c(1,8), col.wdth.ratio=c(4,2),
                            plotAlignments.args = list(line.weight=0.5, ins.size=2, legend.symbol.size=4),
                            plotFreqHeatmap.args = list(plot.text.size=3, x.size=8, legend.text.size=8,
                                                        legend.key.height=grid::unit(0.5, "lines")))
        }
        if (pam == 'Cpf1') {
          p <- plotVariants(crispr_set, txdb=txdb, gene.text.size=8,
                            row.ht.ratio=c(1,8), col.wdth.ratio=c(4,2),
                            plotAlignments.args = list(line.weight=0.5, ins.size=2, legend.symbol.size=4, guide.loc=IRanges(6, 32), pam.start=6, pam.end=9),
                            plotFreqHeatmap.args = list(plot.text.size=3, x.size=8, legend.text.size=8,
                                                        legend.key.height=grid::unit(0.5, "lines")))
        }
        dev.off()

        freqs <- crispr_set$cigar_freqs
        freqs_all <- colSums(freqs)
        indel <- grep('[DI]', rownames(freqs))
        if (length(indel) == 0) {
          freqs_indel <- rep(0, ncol(freqs))
        }
        if (length(indel) == 1) {
          freqs_indel <- freqs[indel,]
        }
        if (length(indel) > 1) {
          freqs_indel <- colSums(freqs[indel,,drop=F])
        }

        indel_percent <- round(freqs_indel / freqs_all * 100, 2)
        sample_names <- names(freqs_all)
        indel_table <- rbind(indel_table, as.data.frame(cbind(plate, gname, sample_names, freqs_all, freqs_indel, indel_percent)))

        snv <- grep('SNV:', rownames(freqs))
        if (length(snv) > 0) {
          snv_freqs <- freqs[snv,,drop=F]
          lst <- strsplit(sub('SNV:', '', rownames(snv_freqs)), ',')
          snv_freqs2 <- NULL
          lst2 <- NULL
          for (n in 1:length(lst)) {
            for (m in 1:length(lst[[n]])) {
              snv_freqs2 <- rbind(snv_freqs2, snv_freqs[n,,drop=F])
              lst2 <- c(lst2, lst[[n]][m])
            }
          }
          snv_freqs3 <- aggregate(snv_freqs2, by=list(lst2), FUN=sum)
          rownames(snv_freqs3) <- snv_freqs3[, 1]
          snv_freqs3 <- snv_freqs3[, -1, drop=F]
          if (ncol(snv_freqs3) == 1) {
            snv_freqs3 <- round(snv_freqs3/freqs_all*100,2)
          } else {
            snv_freqs3 <- t(apply(snv_freqs3, 1, function(x){round(x/freqs_all*100,2)}))
          }
          snv_table[rownames(snv_freqs3), colnames(snv_freqs3)] <- snv_freqs3
        }
        
        DI <- sort(grep('[DI]', rownames(freqs), value = T))
        if (length(DI) > 0){
          DI_freqs <- freqs[DI,,drop=F]
          colnames(DI_freqs) <- sub('._', 'X_', colnames(DI_freqs))
          DI_freqs_all <- data.frame(matrix(nrow=nrow(DI_freqs), ncol=12, data = 0))
          colnames(DI_freqs_all) <- paste0('X_', 1:12)
          DI_freqs_all$RP <- RP
          DI_freqs_all$INDEL <- rownames(DI_freqs)
          DI_freqs_all[, colnames(DI_freqs)] <- DI_freqs
          DI_freqs_all <- DI_freqs_all[, c('RP', 'INDEL', paste0('X_', 1:12))]
          if (nrow(DI_table) == 0){
            DI_table <- DI_freqs_all
          } else {
            DI_table <- rbind(as.matrix(DI_table), as.matrix(DI_freqs_all))
          }
        }

        runs <- sapply(crispr_set$crispr_runs, function(x){length(x$alns)})
        alns[names(runs)] <- runs
      }
      else {
        indel_table <- rbind(indel_table, as.data.frame(cbind(plate, gname, sample_names, freqs_all=0, freqs_indel=0, indel_percent=0)))
      }
    }
    write.table(snv_table, paste0(output_snv, '_', plate, '_', gname, '.tsv'), col.names = NA, quote=F, sep='\t')
    write.table(DI_table, paste0(output_DI, '_', plate, '_', gname, '.tsv'), row.names = F, quote=F, sep='\t')
    
    alns <- data.frame(alns)
    colnames(alns) <- gname
    tot <- cbind(tot, alns)
  }
  write.table(tot, paste0(output_aln, '_', plate, '.tsv'), col.names = NA, quote=F, sep='\t')
}

write.table(indel_table, output_indel, row.names = F, quote=F, sep='\t')

