#####
# Process pair-end reads of barcode-guide-donor cassette to generate a library reference table mapping barcodes to features.
# Cluster barcodes based on sequence similarity and relative abundances.
# 
# OUTPUT FILES
# new_df: data frame containing information on sequence distances for collapsed barcodes
# clusters: low-read barcodes dictionary indicated which barcodes should be clustered
#####

import argparse
import pandas as pd
import editdistance as ed
import numpy as np
import pickle
from collections import Counter
from distance import hamming
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-in', required=True, help="input reference table for barcode clustering", action='store', dest='file')
parser.add_argument('-cutoff', required=True, help="read count cutoff between low and high read barcodes", action='store', dest='cutoff')
parser.add_argument('-segment', required=True, help="current segment to process", action='store', dest='segment')
parser.add_argument('-total_segments', required=True, help="total number of segments to divide into", action='store', dest='total_segments')
args = parser.parse_args()

FILENAME = args.file
ref_table = pd.read_csv(FILENAME)
ref_table['mut_id'] = ["_".join(i.split("_")[:5]) for i in ref_table['sseqid']]

HIGH_LOW_CUTOFF = int(args.cutoff)
high_reads = ref_table[ref_table['read_count'] > HIGH_LOW_CUTOFF]
low_reads = ref_table[ref_table['read_count'] <= HIGH_LOW_CUTOFF]

LENGTH = len(low_reads)
total_segments = int(args.total_segments)
segment = int(args.segment)
start = int((LENGTH/total_segments)*segment)
if (segment+1 == total_segments):
	low_reads = low_reads.iloc[start:]
else:
	stop = int((LENGTH/total_segments)*(segment+1))
	low_reads = low_reads.iloc[start:stop]

BARCODE_CUTOFF = 4
GUIDE_CUTOFF = 2
DONOR_CUTOFF = 8

nl = []

clusters = {}
for i, r in low_reads.iterrows():
	rqseqid = r['qseqid']
	rsseqid = r['mut_id']
	rguide = r['guide']
	rdonor = r['donor']

	candidates = high_reads[(high_reads['mut_id'] == rsseqid)]
	candidates['barcode_distance'] = [hamming(c, rqseqid) for c in candidates['qseqid']]
	if (len(candidates) > 0) and min(candidates['barcode_distance']) <= BARCODE_CUTOFF:
		new_r = candidates[candidates.index==candidates['barcode_distance'].argmin()].iloc[0]
		guide_ed = ed.eval(new_r['guide'], r['guide'])
		donor_ed = ed.eval(new_r['donor'], r['donor'])
		if (guide_ed <= GUIDE_CUTOFF) and (donor_ed <= DONOR_CUTOFF):
			clusters[rqseqid] = new_r['qseqid']
			nl.append([rqseqid, new_r['qseqid'], new_r['barcode_distance'], guide_ed, donor_ed, r['read_count'], new_r['read_count']])

TMP = FILENAME.split("unclustered_reference_table")[0]
new_df = pd.DataFrame(nl, columns=['qseqid1', 'qseqid2', 'barcode_ed', 'guide_ed', 'donor_ed', 'read_count_1', 'read_count_2'])
new_df.to_csv(TMP + str(segment) + "_cluster_details.csv")

pickle_out = open(TMP + str(segment)  + ".cluster_dict", "wb")
pickle.dump(clusters, pickle_out, protocol=2)
pickle_out.close()
print(str(datetime.datetime.now()) + " Finishing segment number " + str(segment))
