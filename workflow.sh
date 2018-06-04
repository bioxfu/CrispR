## Spliting
Rscript script/parse_barcode_in_excel.R
# 
script/split_fastq_by_barcode.py -f raw/JKZ0249/R7-12_T1-1-1/JKZ0249-5583_HK5YMCCXY_L8_1.fq.gz -r raw/JKZ0249/R7-12_T1-1-1/JKZ0249-5583_HK5YMCCXY_L8_2.fq.gz -b tables/JKZ0249.barcode_sequence.tsv -o split/R7-12_T1-1-1
script/split_fastq_by_barcode.py -f raw/JKZ0249/R7-12_T1-1-2/JKZ0249-5584_HK5YMCCXY_L8_1.fq.gz -r raw/JKZ0249/R7-12_T1-1-2/JKZ0249-5584_HK5YMCCXY_L8_2.fq.gz -b tables/JKZ0249.barcode_sequence.tsv -o split/R7-12_T1-1-2
#
script/split_fastq_by_barcode.py -f raw/JKZ0259/R1-R6_T2-1-1/JKZ0259_5671_R1.fq.gz -r raw/JKZ0259/R1-R6_T2-1-1/JKZ0259_5671_R2.fq.gz -b tables/JKZ0259.barcode_sequence.tsv -o split/R1-R6_T2-1-1
script/split_fastq_by_barcode.py -f raw/JKZ0259/R1-R6_T2-1-2/JKZ0259_5672_R1.fq.gz -r raw/JKZ0259/R1-R6_T2-1-2/JKZ0259_5672_R2.fq.gz -b tables/JKZ0259.barcode_sequence.tsv -o split/R1-R6_T2-1-2 

## Mapping
module add bwa/0.7.15
module add samtools/1.3
# build index
cp /home/xfu/Gmatic5/genome/tair10/tair10.fa index/tair10.fa
bwa index index/tair10.fa
# mapping
ls split/ |./script/rush -k 'mkdir -p bam/{}'
find split/*|sed -n 's/_R[12].fastq.gz//p'|sort|uniq|./script/rush -k 'bwa mem index/tair10.fa {}_R1.fastq.gz {}_R2.fastq.gz | samtools view -Shb | samtools sort -o bam/{/%}/{%@split/(.+?)/}.bam'
ls bam/*/*|parallel --gnu 'samtools index {}'

## CripsRVariants
~/R/3.5.0/bin/Rscript script/CripsRVariants.R
