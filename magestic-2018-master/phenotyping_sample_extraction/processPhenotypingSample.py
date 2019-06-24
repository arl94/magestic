import argparse
import gzip
import numpy as np
from collections import Counter
import pandas as pd
from distance import hamming
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f", "-file", help="barcode sequencing fastq file", required=True, action="store", dest="file")
parser.add_argument("-cc", "-cluster_cutoff", required=True, action="store", dest="cluster_cutoff")
parser.add_argument("-bc", "-barcode_cutoff", required=True, action="store", dest="barcode_cutoff")
parser.add_argument("-qc", "-quality_cutoff", default=60, required=False, action="store", dest="quality_cutoff")
parser.add_argument("-o", "-out", required=True, action="store", dest="outfile")
args = parser.parse_args()

# Constants
CLUSTER_CUTOFF = int(args.cluster_cutoff) # cutoff for high/low barcodes during clustering
BARCODE_CUTOFF = int(args.barcode_cutoff) # maximum permitted hamming distance for barcode clustering
QUALITY_CUTOFF = int(args.quality_cutoff) # quality score cutoff for excluding reads
OUTFILE_NAME = str(args.outfile) # Keyword for output files.

# Process barcode sequencing file and filter out low quality reads.
barcode_read_lines = gzip.open(args.file).readlines()
barcode_read_sequence = [barcode_read_lines[r].decode("utf-8").replace("\n","") for r in range(1, len(barcode_read_lines), 4)]
barcode_read_quality = [barcode_read_lines[r].decode("utf-8").replace("\n","") for r in range(3, len(barcode_read_lines), 4)]
barcode_read_quality_scores = [np.mean([ord(i) for i in line]) for line in barcode_read_quality]
barcode_read_sequence = [s for s, qs in zip(barcode_read_sequence, barcode_read_quality_scores) if qs >= QUALITY_CUTOFF]

ref_table = pd.DataFrame.from_dict(Counter(barcode_read_sequence), orient="index") # Make count table.
ref_table.columns = ["read_count"]

# Split barcodes into those with high and low reads for clustering.
high_reads = ref_table[ref_table["read_count"] > CLUSTER_CUTOFF]
low_reads = ref_table[ref_table["read_count"] <= CLUSTER_CUTOFF]
high_reads.sort_values("read_count", ascending=False, inplace=True)

new_barcodes = []
new_barcode_distances = []
high_read_barcodes = list(high_reads.index)

for i, r in low_reads.iterrows(): # Attempt to cluster each low read barcode.
	found = False
	for hb in high_read_barcodes:
		ham = hamming(r.name, hb)
		if ham <= BARCODE_CUTOFF:
			new_barcodes.append(hb)
			new_barcode_distances.append(ham)
			found = True
			break
	if not found:
		new_barcodes.append(r.name)
		new_barcode_distances.append(0)

low_reads["cluster"] = new_barcodes
low_reads["distance"] = new_barcode_distances
high_reads["cluster"] = high_reads.index
high_reads["distance"] = [0 for i in range(0, len(high_reads))]

barcode_df = pd.concat([low_reads, high_reads])
barcode_df.to_csv(OUTFILE_NAME + "_barcode.csv")
cluster_df = barcode_df.groupby(["cluster"])["read_count"].sum()
cluster_df.to_csv(OUTFILE_NAME + "_cluster.csv")
os.system("gzip " + OUTFILE_NAME + "_barcode.csv")
os.system("gzip " + OUTFILE_NAME + "_cluster.csv")
