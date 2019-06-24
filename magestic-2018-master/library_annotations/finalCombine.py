#####
# Process pair-end reads of barcode-guide-donor cassette to generate a library reference table mapping barcodes to features.
# Join all clustered barcode information for individual segments and create a clustered, collapsed library reference table.
# 
# OUTPUT FILES
# ref_table: final, collapsed reference table mapping barcodes to features.
#####

import argparse
import pandas as pd
import numpy as np
import pickle
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-in', required=True, help="mapped cluster dictionaries", nargs='+', action='store', dest='files')
parser.add_argument('-final', required=True, help="combine results file for subpool", action='store', dest='final_file')
args = parser.parse_args()

files = args.files
ff = args.final_file

cluster_dicts = {k: v for d in [pickle.load(open(f, "rb")) for f in files] for k, v in d.items()}
ref_table = pd.read_csv(ff, index_col='qseqid')

cluster = []
for barcode in list(ref_table.index):
	if barcode in cluster_dicts:
		cluster.append(cluster_dicts[barcode])
	else:
		cluster.append(barcode)

ref_table['cluster'] = cluster
labels, levels = pd.factorize(ref_table['cluster'])
ref_table['cluster_id'] = labels
ref_table.to_csv("_".join(ff.split("_")[:2]) + "_final_barcode.csv")
ref_table.sort_values("read_count", ascending=False, inplace=True)
read_counts = pd.DataFrame(ref_table.groupby("cluster_id")['read_count'].sum())
ref_table = ref_table[~ref_table.index.duplicated(keep='first')]
del ref_table['read_count']
ref_table = pd.merge(ref_table, read_counts, left_on='cluster_id', right_index=True, how='left')
ref_table.to_csv("_".join(ff.split("_")[:2]) + "_final_collapsed.csv")


