"""
Process pair-end reads of barcode-guide-donor Step 1 cassette to generate a library reference table mapping barcodes to features.
Map guide-donor consensus sequences to design file using BLAST database, annotate library synthesis errors.
 
OUTPUT FILES
df: Mapped but incomplete reference table for current segment of barcodes.

"""

from Bio import pairwise2
import argparse
import os
import pandas as pd
import gzip
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-in', required=True, help="input filename", action='store', dest='filename')
parser.add_argument('-db', required=True, help="blast database to use", action='store', dest='database')
parser.add_argument('-p', required=True, help="reference plasmid data file (SK1.txt, RM11.txt, etc)", action='store', dest='data_file') #TODO: find a work-around for requiring this argument.
parser.add_argument('-flength', required=False, default=None, help="length of forward read to use in mapping results", action='store', dest='forward_length') #default used to be 140
parser.add_argument('-rlength', required=False, default=None, help="length of reverse read to use in mapping results", action='store', dest='reverse_length')
args = parser.parse_args()

queries = {}
f = open(args.filename).readlines()

OUTPUT_HEADER = args.filename
f = [i.replace("\n", "") for i in f]
for i in range(0, len(f), 2):
	queries[f[i].split("_")[0][1:]] = f[i+1] # ">ATGATG_consensus" to "ATGATG"
#print(str(datetime.datetime.now()) + ' ./blastn -query ' + args.filename + ' -db ' + args.database + ' -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" -out ' + OUTPUT_HEADER + '_blast.csv -num_alignments 1')
df = pd.read_csv(OUTPUT_HEADER.split("consensus")[0] + 'blast.csv', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qseq', 'sseq'])
df = df[['qseqid', 'sseqid', 'pident', 'length']]

data_file = [s.decode('utf-8').replace('\n', '') for s in gzip.open(args.data_file).readlines()]
data = {}
for i in range(0, len(data_file), 2):
	data[data_file[i][1:].replace(",",".")] = data_file[i+1].upper()

ins, dels, mismatches = [], [], []
z1, z2 = [], []

flength = args.forward_length
rlength = args.reverse_length

if flength is not None: flength = int(flength)
else:
	flength = None

if rlength is not None: rlength = int(rlength)

for i, r in df.iterrows():
	curr_ins, curr_dels, curr_mismatches = [], [], []
	qseq = r['qseqid'].split("_")[0]
	sseq = r['sseqid']
	s1 = queries[qseq]
	s2 = data[sseq]
	alignment = pairwise2.align.globalms(s2[:rlength], s1[:flength], 2, -1, -5, -0.5)[0] # flipped order of s2 and s1.
	s1_aln, s2_aln = alignment[1], alignment[0] # flipped order of s2 and s1.
	z1.append(s1_aln)
	z2.append(s2_aln)
	if (s1_aln != s2_aln):
		actual_s2_pos = 0
		new_stop = len(s2_aln)
		for pos in range(0, len(s2_aln)):
			if pos == new_stop: break
			p1 = s1_aln[pos]
			p2 = s2_aln[pos]
			if (p1 != p2) and (p1 != '-') and (p2 != '-'):
				curr_mismatches.append(actual_s2_pos+1)
				actual_s2_pos += 1
			elif (p2 == '-'): # insertion in qseq
				curr_ins.append(actual_s2_pos+1)
			elif (p1 == '-'): # deletion in qseq
				curr_dels.append(actual_s2_pos+1)
				actual_s2_pos += 1
				new_stop -= 1
			else:
				actual_s2_pos += 1 # no mismatch or indel at this position
	ins.append(curr_ins)
	dels.append(curr_dels)
	mismatches.append(curr_mismatches)

df['qseq'] = z1
df['sseq'] = z2
df['insertions'] = ins
df['deletions'] = dels
df['mismatches'] = mismatches
df.to_csv(OUTPUT_HEADER.split("_consensus_queries.fasta")[0] + "_blast_mapped.csv", index=False)
print(str(datetime.datetime.now()) + " Finishing " + OUTPUT_HEADER.split("/")[-1])



