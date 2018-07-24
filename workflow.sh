## Config
PLATE=plate01
SHEET=sample_sheet
READ1=raw/project01/plate01/test_R1.fastq.gz
READ2=raw/project01/plate01/test_R2.fastq.gz
#
BARCODE=tables/$SHEET.barcode_sequence.tsv
BOWTIE_INDEX=~/Gmatic7/genome/tair10/bowtie/tair10
BWA_INDEX=~/Gmatic7/genome/tair10/bwa/tair10
TXDB=~/Gmatic7/gene/tair10/txdb/tair10_txdb.sqlite
RVERSION=3.5.0

## Get barcode information from Excel file
~/R/$RVERSION/bin/Rscript script/parse_barcode_in_excel.R

## Split reads according to the barcode
script/split_fastq_by_barcode.py -f $READ1 -r $READ2 -b $BARCODE -o split/$PLATE

# Map gRNA sequence to genome to find its position
source activate gmatic
bowtie -f -a $BOWTIE_INDEX guide/gRNA.fa |awk '{print $3"\t"$4"\t"$4+length($5)"\t"$1"\t0\t"$2}' > guide/gRNA.bed

# Map short reads
ls split/ |./script/rush -k 'mkdir -p bam/{}'
find split/*|sed -n 's/_R[12].fastq.gz//p'|sort|uniq|./script/rush -k "bwa mem $BWA_INDEX {}_R1.fastq.gz {}_R2.fastq.gz | samtools view -Shb | samtools sort -o bam/{/%}/{%@split/(.+?)/}.bam"
ls bam/*/*|parallel --gnu 'samtools index {}'

## CripsRVariants
~/R/$RVERSION/bin/Rscript script/CripsRVariants.R $TXDB

## Double check the indel frequency
find bam/$PLATE/*.bam -printf "%f\n"|sed 's/.bam//'|parallel --gnu "bedtools bamtobed -i bam/$PLATE/{}.bam -cigar| bedtools intersect -a guide/gRNA.bed -b - -wa -wb |awk '{print \"$PLATE\t\"\$4\"\t{}\t\"\$13}' > {}.tmp"; cat *.tmp > tables/reads_from_gRNA_with_cigar; rm *.tmp
~/R/$RVERSION/bin/Rscript script/indel_from_cigar.R
