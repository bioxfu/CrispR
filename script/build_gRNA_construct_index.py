from Bio import SeqIO
import subprocess

read_len = 150

dict1 = SeqIO.to_dict(SeqIO.parse("guide/construct.fa", "fasta"))
dict2 = SeqIO.to_dict(SeqIO.parse("guide/gRNA.fa", "fasta"))

output = open('guide/gRNA_construct.fa', 'w')
for i,k in enumerate(dict2):
	p5 = str(dict1['5prime'].seq)
	p3 = str(dict1['3prime'].seq)
	gRNA = str(dict2[k].seq)

	p5_len = len(p5)
	p3_len = len(p3)
	gRNA_len = len(gRNA)

	n = read_len - gRNA_len

	if p5_len > n:
		p5_part = p5[-n:]
	else:
		p5_part = p5

	if p3_len > n:
		p3_part = p3[:n]
	else:
		p3_part = p3

	seq = p5_part + gRNA + p3_part
	ID = k + '_' + str(len(p5_part)) + '_' + str(len(gRNA)) + '_' + str(len(p3_part))
	output.write('>' + ID + '\n' + seq + '\n')

output.close()

subprocess.call('bowtie2-build guide/gRNA_construct.fa guide/gRNA_construct.fa', shell=True)
