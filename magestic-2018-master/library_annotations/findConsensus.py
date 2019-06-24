"""
Process pair-end reads of barcode-guide-donor Step 1 cassette to generate a library reference table mapping barcodes to features.
Find consensus guide-donor sequence for each barcode and store in output FASTA file.

results: FASTA file with the consensus sequence for each barcode and corresponding guide-donor sequence.

"""

from Bio import pairwise2
from collections import Counter
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-r1', '-forward', required=True, help="dictionary of barcodes to forward reads", action='store', dest='R1_dict')
parser.add_argument('-r2', '-reverse', required=True, help="dictionary of barcodes to reverse reads", action='store', dest='R2_dict')
parser.add_argument('-o', '-out', required=True, help="keywords to save output files under", action='store', dest='out')
parser.add_argument('-s', '-segment', required=True, help="segment number being processed", action='store', dest='segment') #TODO only used for output file name, should be depreciated
parser.add_argument('-ts', '-total', required=True, help="total number of segments being processed", action='store', dest='total_segments') #TODO only used for output file name, should be depreciated
parser.add_argument('-flength', required=False, default=50, help="length of forward read to use in finding consensus (default=50)", action='store', dest='forward_length')
parser.add_argument('-rlength', required=False, default=None, help="length of reverse read to use in finding consensus (default=None)", action='store', dest='reverse_length')
parser.add_argument('-th', '-threshold', required=False, default=0.50, help="threshold for finding a consensus guide_donor sequence (default=0.50)", action='store', dest='threshold')
args = parser.parse_args()

# Returns the reverse complement DNA sequence.
def reverse_complement(dna):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
	return ''.join([complement[base] for base in dna[::-1]])

# Pads a list of sequences so that they are all the same length.
def create_padded_records(records):
	maxlen = max(len(record) for record in records)
	new_records = []
	for record in records:
		if len(record) != maxlen:
			sequence = str(record).ljust(maxlen, 'N')
			new_records.append(sequence)
		else:
			new_records.append(record)
	return new_records

# Returns a "dumb" consensus sequence for a list of sequences (records)
def get_consensus(new_records, threshold=0.50):
	consensus = ""
	num_records = len(new_records)
	for i in range(0, len(new_records[0])):
		nt = [n[i] for n in new_records]
		count = Counter(nt)
		mc_list = count.most_common()
		most_common = mc_list[0]
		if (len(mc_list) > 1):
			if (mc_list[1][1] == most_common[1]):
				consensus = consensus + "N"
				continue
		mc_nt = most_common[0]
		mc_pct = most_common[1] / float(num_records)
		if mc_pct >= threshold:
			consensus = consensus + mc_nt
		else:
			consensus = consensus + "N"
	return consensus

# Generates the merged consensus sequence for a particular barcode.
def merge_consensus(consensus2, short_consensus1):
	alignment = pairwise2.align.globalms(consensus2, short_consensus1, 2, -1, -5, -0.5)[0]
	first_aln = alignment[0]
	second_aln = alignment[1]

	merged = ""
	for i in range(0, len(first_aln)):
		first = first_aln[i]
		second = second_aln[i]
		if first != "-":
			merged = merged + first
		elif second != "-":
			merged = merged + second
	return merged

OUTPUT_HEADER = args.out
THRESHOLD = float(args.threshold)

segment = int(args.segment)
total_segments = int(args.total_segments)

pickle_in = open(args.R1_dict, "rb")
R1_dict = pickle.load(pickle_in)
pickle_in = open(args.R2_dict, "rb")
R2_dict = pickle.load(pickle_in)

barcodes = list(R1_dict.keys())
results = open(OUTPUT_HEADER + "_" + str(segment) + "-" + str(total_segments) + '_consensus_queries.fasta', 'w')

flength = args.forward_length
rlength = args.reverse_length

for b in range(0, len(barcodes)):
	bar = barcodes[b]
	records1 = R1_dict[bar]
	records2 = R2_dict[bar]

	new_records1 = create_padded_records(records1)
	new_records2 = create_padded_records(records2)
	consensus1 = reverse_complement(get_consensus(new_records1, threshold=THRESHOLD))
	consensus2 = get_consensus(new_records2, threshold=THRESHOLD)
	
	if flength is not None:
		consensus1 = consensus1[:int(flength)]
	if rlength is not None:
		consensus2 = consensus2[:int(rlength)]

	merged = merge_consensus(consensus2, consensus1)
	results.write(">" + bar + "_consensus\n")
	results.write(merged + "\n")
results.close()