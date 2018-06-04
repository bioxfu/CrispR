library(GenomicFeatures)

gtf_fname <- './txdb/tair10.gtf'

txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = 'gtf')

saveDb(txdb, file='./txdb/tair10_txdb.sqlite')

