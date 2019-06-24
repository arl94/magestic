import numpy as np
import argparse
import pandas as pd
import analysis

parser = argparse.ArgumentParser()
parser.add_argument('-t', required=True, help="prepared counts table file", action='store', dest='counts_table')
parser.add_argument('-m', required=True, help="model designed", action='store', dest='model_design')
parser.add_argument('--barcode', required=False, const=True, nargs="?", help="barcode level analysis indicator", action='store', dest='barcode_level')
parser.add_argument('--feature', required=False, const=True, nargs="?", help="feature level analysis indicator", action='store', dest='feature_level')
parser.add_argument('--v2', required=False, const=True, nargs="?", help="v1 feature analysis indicator", action='store', dest='v2')
parser.add_argument('-o', required=True, help="desired path/name of output analysis file", action='store', dest='out')
args = parser.parse_args()

counts = pd.read_csv(args.counts_table, index_col=0)
model = pd.read_csv(args.model_design, index_col=0)
#model = model.to_numeric()
cols = model.columns
for col in cols:
	model[col] = model[col].astype(int)

outfile = args.out.split(".csv")[0]

samples = list(model.index)
#counts['library'] = counts['library'].fillna("unknown")
#counts['subpool'] = counts['subpool'].fillna("unknown")
counts['guide_status'] = counts['guide_status'].fillna("unknown")
#counts['donor_status'] = counts['donor_status'].fillna("unknown")
print(samples)
if args.barcode_level:
	counts['unique_id'] = list(counts.index)
	label = "barcode_level"
elif args.feature_level:
	counts['unique_id'] = [str(i)+"_"+str(j)+"_"+str(k)+"_"+str(l) for i,j,k,l in zip(counts['sseqid'], counts['guide_status'], counts['donor_status'], counts['library'])]
	label = "feature_level"
"""elif args.v2:
	counts['unique_id'] = [str(i)+"_"+str(j)+"_"+str(k) for i,j,k in zip(counts['ref'], counts['guide_status'], counts['MD'])]
	label = "feature_level"	
else: # default is barcode level analysis
	counts['unique_id'] = list(counts.index)
	label = "barcode_level"
"""
#counts = counts[samples + ['unique_id', 'ref', 'guide_status', 'MD']].dropna(subset=samples, how='all')
counts = counts[samples + ["unique_id", "sseqid", "library", "subpool", "guide_status", "donor_status", "bsp_status"]].dropna(subset=samples, how='all')
counts['percentNA'] = counts[samples].isnull().sum(axis=1) / float(len(samples))
counts = counts[counts['percentNA'] < 0.25]

counts.to_csv(outfile + "_" + label + "_counts.csv")
results, log = analysis.run(counts, model, log=True)
results.to_csv(outfile + "_" + label + "_analysis.csv")
log.to_csv(outfile + "_" + label + "_log.csv")
