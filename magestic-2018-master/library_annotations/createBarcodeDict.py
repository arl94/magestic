"""
Process pair-end reads of barcode-guide-donor Step 1 cassette to generate a library reference table mapping barcodes to features.
Create dictionaries mapping barcodes to forward and reverse reads, split into sub-segments.

R1_dict: map barcodes to corresponding R1 sequences.
R2_dict: map barcodes to corresponding R2 sequences.
read_count_dict: map each barcode to corresponding total number of reads.

"""

from collections import Counter
import argparse
import gzip
import numpy as np
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-f', '-forward', required=True, help="forward sequencing files", nargs='+', action='store', dest='forward_files')
parser.add_argument('-r', '-reverse', required=True, help="reverse sequencing files", nargs='+', action='store', dest='reverse_files')
parser.add_argument('-s', '-segments', required=True, help="number of segments to split job into", action='store', dest='total_segments')
parser.add_argument('-o', '-out', required=True, help="keyword for saving output files", action='store', dest='out')
parser.add_argument('-c', '-cutoff', required=False, default=0, help="read count cutoff for barcodes to keep (default=0)", action='store', dest='cutoff')
parser.add_argument('-b', '-barcode', required=False, default=31, help="length of barcode (default=31)", action='store', dest='barcode_length')
parser.add_argument('-bq', '-bquality', required=False, default=53, help="ascii quality score cutoff for barcode (default=53)", action='store', dest='barcode_quality')
parser.add_argument('-gdq', '-gdquality', required=False, default=55, help="ascii quality score cutoff for guide-donor (default=55)", action='store', dest='guide_donor_quality')

args = parser.parse_args()

OUTPUT_HEADER = args.out
READ_COUNT_CUTOFF = int(args.cutoff)
BARCODE_LENGTH = int(args.barcode_length)
BARCODE_QUALITY_CUTOFF = int(args.barcode_quality)
GUIDE_DONOR_QUALITY_CUTOFF = int(args.guide_donor_quality)

# Collect all sequencing reads from forward files.
forward_lines = []
for file in args.forward_files:
	forward_lines.extend(gzip.open(file).readlines())

# Forward sequence.
forward_sequence = [forward_lines[r] for r in range(1, len(forward_lines), 4)]
forward_sequence = [l.decode('utf-8').replace("\n","") for l in forward_sequence]

# Forward sequence quality scores.
forward_quality = [forward_lines[r] for r in range(3, len(forward_lines), 4)]
forward_quality = [l.decode('utf-8').replace("\n","") for l in forward_quality]

barcode_quality_scores = [] # Barcode quality.
for line in forward_quality:
        scores = [ord(i) for i in line[:BARCODE_LENGTH]]
        barcode_quality_scores.append(np.mean(scores))

forward_guide_donor_quality_scores = [] # Guide-donor quality.
for line in forward_quality:
        scores = [ord(i) for i in line[BARCODE_LENGTH:]]
        forward_guide_donor_quality_scores.append(np.mean(scores))

# Collect all sequencing reads from reverse files.
reverse_lines = []
for file in args.reverse_files:
	reverse_lines.extend(gzip.open(file).readlines())

# Reverse sequence.
reverse_sequence = [reverse_lines[r] for r in range(1, len(reverse_lines), 4)]
reverse_sequence = [l.decode('utf-8').replace("\n","") for l in reverse_sequence]

# Reverse sequence base quality scores.
reverse_quality = [reverse_lines[r] for r in range(3, len(reverse_lines), 4)]
reverse_quality = [l.decode('utf-8').replace("\n","") for l in reverse_quality]

reverse_guide_donor_quality_scores = []
for line in reverse_quality:
        scores = [ord(i) for i in line]
        reverse_guide_donor_quality_scores.append(np.mean(scores))

# Filter out low quality barcodes and low quality guide-donor sequences.
forward_sequence, reverse_sequence, barcodes = zip(*[(f, r, f[:BARCODE_LENGTH]) for f, r, fscore, fscore2, rscore
        in zip(forward_sequence, reverse_sequence, barcode_quality_scores, forward_guide_donor_quality_scores, reverse_guide_donor_quality_scores) 
        if (fscore >= BARCODE_QUALITY_CUTOFF) and (fscore2 >= GUIDE_DONOR_QUALITY_CUTOFF) and (rscore >= GUIDE_DONOR_QUALITY_CUTOFF)])

if (READ_COUNT_CUTOFF != 0): # optional choice to remove low read barcodes from annotations.
	barcodes_to_keep = [key for key, count in Counter(barcodes).items() if count >= READ_COUNT_CUTOFF]
	keep_dict = {g: True for g in barcodes_to_keep}
	forward_sequence, reverse_sequence, barcodes = zip(*[(f, r, b) for f, r, b 
		in zip(forward_sequence, reverse_sequence, barcodes) if b in keep_dict])

# Store barcode read count dictionary for later use. 
count_dict = dict(Counter(barcodes))
pickle_out = open(OUTPUT_HEADER + ".read_count_dict", "wb")
pickle.dump(count_dict, pickle_out, protocol=2)
pickle_out.close()

# Divide up barcodes into specified number of segments for parallel analysis.
LENGTH = len(set(barcodes))
total_segments = int(args.total_segments)

barcode_list = list(set(barcodes))
for segment in range(0, total_segments):
	start = int((LENGTH/total_segments)*segment) # determine start and end position of segment.
	if (segment+1 == total_segments):
		sub_barcodes_set = barcode_list[start:]
	else:
		stop = int((LENGTH/total_segments)*(segment+1))
		sub_barcodes_set = barcode_list[start:stop]
	sub_barcodes_dict = {b: True for b in sub_barcodes_set}

	sub_forward, sub_reverse, sub_barcodes = zip(*[(f, r, b) for f, r, b 
		in zip(forward_sequence, reverse_sequence, barcodes) if b in sub_barcodes_dict])

	R1_dict, R2_dict = {}, {} # store reads by barcode into R1 and R2 dictionaries.
	for f, r, b in zip(sub_forward, sub_reverse, sub_barcodes):
		if (b not in R1_dict) and (b not in R2_dict):
			R1_dict[b] = [f]
			R2_dict[b] = [r]
		else:
			R1_dict[b].append(f)
			R2_dict[b].append(r)

	pickle_out = open(OUTPUT_HEADER + "_" + str(segment) + "-" + str(total_segments) + ".R1_dict", "wb")
	pickle.dump(R1_dict, pickle_out, protocol=2)
	pickle_out.close()

	pickle_out = open(OUTPUT_HEADER + "_" + str(segment) + "-" + str(total_segments) + ".R2_dict", "wb")
	pickle.dump(R2_dict, pickle_out, protocol=2)
	pickle_out.close()