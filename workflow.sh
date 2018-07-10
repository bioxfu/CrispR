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

# gRNA position
module add bowtie/1.1.2
bowtie -f -v 0 -a ~/Gmatic5/genome/tair10/tair10 guide/gRNA_R1-R6.fa |awk '{print $3"\t"$4"\t"$4+length($5)"\t"$1"\t0\t"$2}' > guide/gRNA_R1-R6.bed
bowtie -f -v 0 -a ~/Gmatic5/genome/tair10/tair10 guide/gRNA_R19-R24.fa |awk '{print $3"\t"$4"\t"$4+length($5)"\t"$1"\t0\t"$2}' > guide/gRNA_R19-R24.bed

## CripsRVariants
~/R/3.5.0/bin/Rscript script/CripsRVariants.R

