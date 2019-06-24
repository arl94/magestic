"""
Process pair-end reads of barcode-guide-donor Step 1 cassette to generate a library reference table mapping barcodes to features.
Join all reference table information for individual segments, add read count information and guide/donor status annotations.

df: uncollapsed reference table mapping barcodes to features.

"""

import pickle
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', required=True, help="mapped segment files", nargs='+', action='store', dest='files')
parser.add_argument('-c', required=True, help="count dict", action='store', dest='count_file')
parser.add_argument('-l', required=True, help="name of library", action='store', dest='library')
parser.add_argument('-s', required=True, help="number of segments to combine", action='store', dest='segs')
args = parser.parse_args()

files = args.files
count_file = args.count_file
SEGS = int(args.segs)

frames = []
for i in range(0, SEGS):
	frames.append(pd.read_csv(files[i]))
result = pd.concat(frames)

pickle_in = open(count_file, "rb")
count_dict = pickle.load(pickle_in)
count_df = pd.DataFrame.from_dict(count_dict, orient='index')
count_df.columns = ['read_count']

result['qseqid'] = [i.split("_")[0] for i in list(result['qseqid'])]
result = pd.merge(result, count_df, left_on="qseqid", right_index=True, how="left")
subpool = "_".join([files[0].split("_")[0], files[0].split("_")[1]])
result['subpool'] = subpool.split("_")[0]
result['guide'] = [i[:20] for i in result['qseq']]
result['donor'] = [i[31:] for i in result['qseq']]

result2 = result[(result['length'] >= 125) & (result['length'] <= 155)]

#TODO: add T-score annotations!

# Annotate guide status.
guide_status_list = []
donor_status_list = []
bsp_status_list = []
for i, r in result.iterrows():
	indels = eval(r['insertions']) + eval(r['deletions'])
	mismatches = eval(r['mismatches'])

	guide_indels, guide_mismatches = [], []
	donor_indels, donor_mismatches = [], []
	bsp_indels, bsp_mismatches = [], []
	for i in indels:
		if i <= 20:
			guide_indels.append(i)
		elif 20 < i <= 31:
			bsp_indels.append(i)
		else: # i > 31
			donor_indels.append(i)

	for m in mismatches:
		if m <= 20:
			guide_mismatches.append(m)
		elif 20 < m <= 31:
			bsp_mismatches.append(m)
		else: # m > 31
			donor_mismatches.append(m)

	curr_dead_guide = 0
	curr_near_perfect_guide = 0
	curr_intermediate_guide = 0

	if (len(guide_indels) != 0):
		for ind in guide_indels:
			if 20 >= int(ind) >= 6:
				curr_dead_guide = 1
			elif 5 >= int(ind) >= 3:
				curr_intermediate_guide = 1
			elif 2 >= int(ind) >= 1:
				curr_near_perfect_guide = 1
	mismatches_in_1_to_15 = 0
	if (len(guide_mismatches) != 0):
		for mm in guide_mismatches:
			if 20 >= int(mm) >= 6:
				mismatches_in_1_to_15 += 1
			elif 5 >= int(mm) >= 3:
				curr_intermediate_guide = 1
			elif (2 >= int(mm) >= 1) and (not curr_intermediate_guide):
				curr_near_perfect_guide = 1
		if mismatches_in_1_to_15 >= 2:
			curr_dead_guide = 1

	if curr_dead_guide == 1:
		guide_status_list.append("dead_guide")
	elif curr_intermediate_guide == 1:
		guide_status_list.append("intermediate_guide")
	elif curr_near_perfect_guide == 1:
		guide_status_list.append("near_perfect_guide")
	elif (len(donor_mismatches) == 0) and (len(donor_indels) == 0) and (len(bsp_mismatches) == 0) and (len(bsp_indels) == 0):
		guide_status_list.append("perfect_guide_donor")
	elif (len(donor_mismatches) == 0) and (len(donor_indels) == 0):
		guide_status_list.append("perfect_guide_donor_imperfect_bsp")
	else:
		guide_status_list.append("perfect_guide")

	if (len(donor_mismatches) == 0) and (len(donor_indels) == 0):
		donor_status_list.append("perfect_donor")
	else:
		donor_status_list.append("imperfect_donor")

	if (len(bsp_mismatches) == 0) and (len(bsp_indels) == 0):
		bsp_status_list.append("perfect_bspQI")
	else:
		bsp_status_list.append("imperfect_bspQI")

result['guide_status'] = guide_status_list
result['donor_status'] = donor_status_list
result['bsp_status'] = bsp_status_list
result['library'] = args.library
result.to_csv(subpool + "_final.csv", index=False)
