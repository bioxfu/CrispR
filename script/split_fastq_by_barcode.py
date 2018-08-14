#! ~/miniconda2/bin/python

## Note: This script only works in Python2

from optparse import OptionParser
from itertools import izip
import gzip
import os

parser = OptionParser()
parser.add_option('-f', '--read1')
parser.add_option('-r', '--read2')
parser.add_option('-b', '--barcode', help='run parse_barcode_in_excel.R to get the barcode file')
parser.add_option('-o', '--output', help='output directory')

(options, args) = parser.parse_args()
if options.read1 is None or options.read2 is None or options.barcode is None or options.output is None:
	parser.error("incorrect number of arguments")  

if not os.path.isdir(options.output):
	os.makedirs(options.output)

fileA = gzip.open(options.read1)
fileB = gzip.open(options.read2)

barcode2sample = {}
with open(options.barcode) as f:
	for line in f:
		RP, FP, RP_seq, FP_seq = line.strip().split('\t')
		if RP_seq != 'RP_sequence':
			barcode2sample[RP_seq+'_'+FP_seq] = RP+'_'+FP

barcode2counts = {}

# open output file
sample2handle = {}
for samples in barcode2sample.values():
	sample2handle[samples] = [gzip.open(options.output+'/'+samples+'_R1.fastq.gz', 'wb'), gzip.open(options.output+'/'+samples+'_R2.fastq.gz', 'wb')]

for lineA, lineB in izip(fileA, fileB):
	lineA_1 = lineA.strip()
	lineA_2 = fileA.next().strip()
	lineA_3 = fileA.next().strip()
	lineA_4 = fileA.next().strip()
	lineB_1 = lineB.strip()
	lineB_2 = fileB.next().strip()
	lineB_3 = fileB.next().strip()
	lineB_4 = fileB.next().strip()

	FP = lineA_2[4:8]
	RP = lineB_2[4:8]
	barcode = RP+'_'+FP
	if barcode in barcode2sample:
		sample = barcode2sample[barcode]
		sample2handle[sample][0].write("%s\n%s\n%s\n%s\n" % (lineA_1+' '+FP, lineA_2[27:], lineA_3, lineA_4[27:]))
		sample2handle[sample][1].write("%s\n%s\n%s\n%s\n" % (lineB_1+' '+RP, lineB_2[27:], lineB_3, lineB_4[27:]))
	
		if barcode in barcode2counts:
			barcode2counts[barcode] += 2
		else:
			barcode2counts[barcode] = 2

# close output file
for samples in barcode2sample.values():
	sample2handle[samples][0].close()
	sample2handle[samples][1].close()

with open(options.output+'/reads_stat.tsv', 'w') as f:
	for barcode in sorted(barcode2counts.keys()):
		f.write("%s\t%s\t%s\n" % (barcode2sample[barcode], barcode, barcode2counts[barcode]))
