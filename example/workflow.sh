## User config
SHEET=sample_sheet # without file extension name (.xlsx)
GRNA=guide/gRNA_ath # without file extension name (.fa)
PLATES=(5686)
OUTPUT_INDEL=tables/indel_freq_all_ath.tsv
OUTPUT_SNV=tables/snv_freq_all_ath # without any file extension name
OUTPUT_ALN=tables/aln_freq_all_ath # without any file extension name

## Other config
MYHOME=/cluster/home/xfu
RVERSION=3.5.1

## Arabidopsis thaliana
FASTA=$MYHOME/Gmatic7/genome/tair10/tair10.fa
BOWTIE_INDEX=$MYHOME/Gmatic7/genome/tair10/bowtie/tair10
BWA_INDEX=$MYHOME/Gmatic7/genome/tair10/bwa/tair10
TXDB=$MYHOME/Gmatic7/gene/tair10/txdb/tair10_txdb.sqlite

## Solanum lycopersicum
# FASTA=$MYHOME/Gmatic7/genome/tomato/Sly3.fa
# BOWTIE_INDEX=$MYHOME/Gmatic7/genome/tomato/bowtie/Sly3
# BWA_INDEX=$MYHOME/Gmatic7/genome/tomato/bwa/Sly3
# TXDB=$MYHOME/Gmatic7/gene/tomato/txdb/Sly3_txdb.sqlite

## Morus notabilis
# FASTA=$MYHOME/Gmatic7/genome/morus/index/morus_notabilis.fa
# BOWTIE_INDEX=$MYHOME/Gmatic7/genome/morus/index/morus_notabilis
# BWA_INDEX=$MYHOME/Gmatic7/genome/morus/index/morus_notabilis
# TXDB=$MYHOME/Gmatic7/genome/morus/txdb/morus_notabilis_txdb.sqlite

## Workflow
echo 'Get barcode information from Excel file'
$MYHOME/R/$RVERSION/bin/Rscript script/parse_barcode_in_excel.R

echo 'Map gRNA sequence to genome to find its position'
if [ ! -f ${GRNA}_fix.bed ]; then
	bowtie -f -v 0 -a $BOWTIE_INDEX ${GRNA}.fa |awk '{print $3"\t"$4"\t"$4+length($5)"\t"$1"\t0\t"$2}' > ${GRNA}.bed
	## sometimes the length of given gRNA is not 23nt, we need to fix the bed file
	grep '+$' ${GRNA}.bed|awk 'BEGIN { OFS = "\t" } {$2=$3-23; print $0}' > ${GRNA}_fix.bed
	grep '\-$' ${GRNA}.bed|awk 'BEGIN { OFS = "\t" } {$3=$2+23; print $0}' >> ${GRNA}_fix.bed
fi

echo 'Split reads according to the barcode'
for PLATE in ${PLATES[@]}; do
	echo $PLATE
	if [ ! -f split/$PLATE/reads_stat.tsv ]; then
	    $MYHOME/miniconda2/bin/python script/split_fastq_by_barcode.py -f fastq/${PLATE}_R1.fastq.gz -r fastq/${PLATE}_R2.fastq.gz -b tables/${SHEET}.barcode_sequence.tsv -o split/$PLATE
	fi
done

echo 'Map short reads'
for PLATE in ${PLATES[@]}; do
	echo $PLATE
	mkdir -p bam/${PLATE}
	find split/${PLATE}/*|sed -n 's/_R[12].fastq.gz//p'|sort|uniq|./script/rush -k "bwa mem $BWA_INDEX {}_R1.fastq.gz {}_R2.fastq.gz | samtools view -Shb | samtools sort -o bam/{/%}/{%@split/(.+?)/}.bam"
	ls bam/${PLATE}/*.bam|parallel --gnu 'samtools index {}'
done

echo 'CripsRVariants'
rm -f ${GRNA}.tsv
for PLATE in ${PLATES[@]}; do  echo -e "$PLATE\t${GRNA}_fix.bed" >> ${GRNA}.tsv; done
$MYHOME/R/$RVERSION/bin/Rscript script/CripsRVariants.R $TXDB $FASTA ${GRNA}.tsv $OUTPUT_INDEL $OUTPUT_SNV $OUTPUT_ALN

#echo 'Double check the indel frequency'
#find bam/$PLATE/*.bam -printf "%f\n"|sed 's/.bam//'|parallel --gnu "bedtools bamtobed -i bam/$PLATE/{}.bam -cigar| bedtools intersect -a guide/gRNA.bed -b - -wa -wb |awk '{print \"$PLATE\t\"\$4\"\t{}\t\"\$13}' > {}.tmp"; cat *.tmp > tables/reads_from_gRNA_with_cigar; rm *.tmp
#$MYHOME/R/$RVERSION/bin/Rscript script/indel_from_cigar.R
