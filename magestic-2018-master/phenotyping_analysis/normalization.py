import pandas as pd
import numpy as np

def __meanColumnSum(df):
	col_sums = list(df.sum())
	mean_col_sum= np.mean(col_sums)
	col_factors = col_sums/mean_col_sum
	return df/col_factors

def __modeNormalization(df, log=False):
	ref_table = np.log2(df)
	means = list((np.mean(ref_table.T)))
	avg_ratios = (ref_table.T-means).T # a - rowavg(a+...+e')
	column_modes = [__half_sample_mode(avg_ratios[c].dropna()) for c in avg_ratios.columns]
	ref_table = 2 ** (ref_table - column_modes)
	return ref_table, column_modes

# Returns estimated mode of continous data distribution.
def __half_sample_mode(x, already_sorted=False):
    if len(x) < 3:
        return np.mean(x)
    if already_sorted:
        sorted_x = x # No need to sort
    else:
        sorted_x = np.sort(x)
    half_idx = int((len(sorted_x) + 1) / 2) # Round up to include the middle value, in the case of an odd-length array

    # Calculate all interesting ranges that span half of all data points
    ranges = sorted_x[-half_idx:] - sorted_x[:half_idx]
    smallest_range_idx = np.argmin(ranges)

    # Now repeat the procedure on the half that spans the smallest range
    x_subset = sorted_x[smallest_range_idx : (smallest_range_idx+half_idx)]
    return __half_sample_mode(x_subset, already_sorted=True)

# Split count table into count-only table and annotations-only table.
def preprocessDefaultExperiment(table_name, samples, col, sep_norm=False): #TODO sampels should be a model_design df..take index values for sample names
	ref_table = pd.read_csv(table_name, index_col=col)
	guide_annotations = ref_table[["sseqid", "library", "subpool", "guide_status", "donor_status", "bsp_status"]]
	if sep_norm:
		rm11_counts = ref_table[ref_table['library']=='RM11_V3'][samples] 
		sk1_counts = ref_table[ref_table['library']=='SK1_V3'][samples] 
		return rm11_counts, sk1_counts, guide_annotations
	counts = ref_table[samples]
	return counts, guide_annotations

#  Remove any rows where the row mean is less than the threshold.
def removeLowAbundantCounts(ref_table, threshold=20):
	return ref_table[~(np.mean(ref_table.T) < threshold)] 

# Remove barcodes that don't have designed TG sites.
def sanitize_barcodes(ref_table, returnDirtyBarcodes=False):
	barcodes = list(ref_table.index)
	good_barcodes, bad_barcodes = [], []
	for b in barcodes:
		if (((b[3:5] == "TG") & (b[10:12] == "TG") & (b[17:19] == "TG") & (b[24:26] == "TG")) or 
		((b[5:7] == "TG") & (b[12:14] == "TG") & (b[19:21] == "TG") & (b[26:28] == "TG"))):
			good_barcodes.append(b)
		else:
			bad_barcodes.append(b)
	if returnDirtyBarcodes:
		bad_barcode_table = ref_table[ref_table.index.isin(bad_barcodes)]
		ref_table = ref_table[ref_table.index.isin(good_barcodes)]
		return ref_table, bad_barcode_table
	else:
		ref_table = ref_table[ref_table.index.isin(good_barcodes)]
		return ref_table

def normalize(df, method="meanColumnSum"):
	dispatcher = {"meanColumnSum": __meanColumnSum, "modeNormalization": __modeNormalization}
	return dispatcher[method](df)






def filter(ref_table, model_design):
	t0 = model_design[model_design['time'] == 0]
	mm1 = model_design[(model_design['time'] == 1) and (model_design['perturbation'] == 0)]
	mm2 = model_design[(model_design['time'] == 2) and (model_design['perturbation'] == 0)]
	x1 = model_design[(model_design['time'] == 1) and (model_design['perturbation'] == 1)]
	x2 = model_design[(model_design['time'] == 1) and (model_design['perturbation'] == 1)]

	cols = t0 + mm1 + mm2 + x1 + x2 #+ y1 + y2
	cols = [c for c in cols if c != "None"]
	ref_table = ref_table[cols] # Select only specified columns for normalization.
	ref_table = ref_table.dropna(subset=[cols], how='all') # drop rows that have all missing values.
	
	ref_table = sanitize_barcodes(ref_table)

	#bmeans = np.mean(bad_barcode_table[cols].T)
	#bmeans = np.log2(bmeans)
	#bnulls = bad_barcode_table[cols].isnull().sum(axis=1)

	ref_table = ref_table.replace(0, np.nan) # Fill all zeros with NAs.
	ref_table_unfiltered = ref_table.copy()
	ref_table_unfiltered = np.log2(ref_table_unfiltered[cols])
	
	print("Processing missing values...")
	before = len(ref_table)
	if (mm1[0] != 'None'):
		ref_table = ref_table[((ref_table[mm1].isnull().sum(axis=1)) <= 1)] # Max. 1/5 NAs for experimental time point without any pertubation.
	if (mm2[0] != 'None'):
		ref_table = ref_table[((ref_table[mm2].isnull().sum(axis=1)) <= 1)] # Max. 1/5 NAs for experimental time point without any pertubation.

	# Exclude rows with 1+ NAs in pertubation and <50 counts in 1+ reps at experimental time point without any pertubation.
	if (x1[0] != 'None'):
		ref_table = ref_table[~(((ref_table[x1].isnull().sum(axis=1)) > 2)  
			& (((ref_table[mm1]).T.replace(np.nan, 50).min()) < 50))]

	if (x2[0] != 'None'):
		ref_table = ref_table[~(((ref_table[x2].isnull().sum(axis=1)) > 2)  
			& (((ref_table[mm2]).T.replace(np.nan, 50).min()) < 50))]

	'''if (y1[0] != 'None'):
		ref_table = ref_table[~(((ref_table[y1].isnull().sum(axis=1)) > 1)  
			& (((ref_table[mm1]).T.replace(np.nan, 50).min()) < 50))]

	if (y2[0] != 'None'):
		ref_table = ref_table[~(((ref_table[y2].isnull().sum(axis=1)) > 1)  
			& (((ref_table[mm2]).T.replace(np.nan, 50).min()) < 50))]
	'''
	print(str(before - len(ref_table)) + " rows removed.")




def preprocessV2Experiment(table_name):
	ref_table = pd.read_csv(table_name, index_col='qseqid')
	del ref_table['Unnamed: 0']
	guide_annotations = ref_table[["MD", "sseqid", "guide_status"]]
	ref_table = ref_table.loc[:, ref_table.columns != 'MD']
	ref_table = ref_table.loc[:, ref_table.columns != 'sseqid']
	ref_table = ref_table.loc[:, ref_table.columns != 'guide_status']
	ref_table = ref_table.fillna(0)
	ref_table = nz.removeLowAbundantCounts(ref_table)
	ref_table = ref_table.replace(0, np.nan)
	ref_table, modes = nz.normalize(ref_table, "modeNormalization")

	ref_table = pd.merge(ref_table, guide_annotations, left_index=True, right_index=True, how='left') # Merge files on barcodes.
	ref_table = ref_table[~ref_table.index.duplicated(keep=False)]

	unique_ids = []
	for i, r in ref_table.iterrows():
		if (r['guide_status'] in ['dead_guide', 'intermediate_guide', 'near_perfect_guide', 'unknown']):
			ref = r['sseqid'] +"_"+r['guide_status']+"_"+i
		else:
			ref = r['sseqid'] +"_"+r['guide_status']
		unique_ids.append(ref)
	ref_table['unique_id'] = unique_ids
	ref_table.sort_values(by=['unique_id'], inplace=True) # Sort table by reference annotation.









"""
table_name = args.counts_table
print("Loading counts table file..." + str(table_name))

# Check whether files are in csv or txt format.
ref_table = loadFile(table_name, ic=0)
guide_annotations = ref_table[["sseqid", "library", "subpool", "guide_status"]]
t0, mm1, mm2, x1, x2, y1, y2 = args.t0, args.mm1, args.mm2, args.x1, args.x2, args.y1, args.y2
cols = t0 + mm1 + mm2 + x1 + x2 + y1 + y2
cols = [c for c in cols if c != "None"]

ref_table = ref_table[cols + ['library']]

# Optional shuffling
if (args.shuffle):
	print("Shuffling data columns...")
	ref_table.columns = random.sample(list(ref_table.columns[:-1]), len(ref_table.columns)-1) + ['library']
	ref_table = ref_table.reindex_axis(sorted(ref_table.columns), axis=1)

# Split dataset into library-specific subsets.
rm11 = ref_table[ref_table['library']== 'RM11_V3']
sk1 = ref_table[ref_table['library']== 'SK1_V3']

# Mode-based normalization.
rm11_norm, rm11_modes = normalization(rm11, t0, mm1, mm2, x1, x2, y1, y2, 'RM11_V3')
sk1_norm, sk1_modes = normalization(sk1, t0, mm1, mm2, x1, x2, y1, y2, 'SK1_V3')

modes = pd.DataFrame(columns=cols)
modes.loc[len(modes)] = rm11_modes
modes.loc[len(modes)] = sk1_modes
modes.index = ["RM11_V2", "SK1_V3"]
modes.to_csv(args.out+"_modes.csv")
ref_table = pd.concat([rm11_norm, sk1_norm])

before = len(ref_table)
ref_table = pd.merge(ref_table, guide_annotations, left_index=True, right_index=True, how='left') # Merge files on barcodes.
ref_table = ref_table[~ref_table.index.duplicated(keep=False)]

print("Creating unique identifiers...")
unique_ids = []
for i, r in ref_table.iterrows():
	if (r['guide_status'] in ['dead_guide', 'intermediate_guide', 'near_perfect_guide', 'unknown']):
		ref = r['sseqid'] +"_"+r['library']+"_"+r['guide_status']+"_"+i
	else:
		ref = r['sseqid'] +"_"+r['library']+"_"+r['guide_status']
	unique_ids.append(ref)
ref_table['unique_id'] = unique_ids
ref_table.sort_values(by=['unique_id'], inplace=True) # Sort table by reference annotation.

print("Generating information on misssing values...")
ref_table['rowmean'] = ref_table[cols].T.mean()
ref_table['rownulls'] = ref_table[cols].isnull().sum(axis=1)

r = str(len(ref_table[~pd.isnull(ref_table['sseqid'])]))
print("Number of barcodes: " + str(len(ref_table)))
print("Number of annotated barcodes: " + r)
print("Number of unique features: " + str(len(ref_table['sseqid'].value_counts())))
print("Number of unique identifiers: " + str(len(ref_table['unique_id'].value_counts())))
ref_table.to_csv(args.out) # save count table to a CSV for analysis.
"""








