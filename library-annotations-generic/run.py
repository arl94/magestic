import argparse
import os
import gzip
import pickle
import numpy as np
from collections import Counter
import multiprocessing
import subprocess
import pandas as pd
import datetime

parser = argparse.ArgumentParser()

# Required command line arguments
parser.add_argument("-d", "-dir",
	required = True,
	help = "name of library design directory",
	action = "store",
	dest = "path")

parser.add_argument("-r", "-ref",
	required=True,
	help="reference guide-donor design file",
	action='store',
	dest='design_file')

parser.add_argument('-db', "-database",
	required=True,
	help="blastn database to use for mapping",
	action='store',
	dest='database')

# Optional arguments
parser.add_argument("-c", "-cpu",
	required=False,
	default=int((multiprocessing.cpu_count()/2.0)),
	help="the number of partitions to split data into (i.e. num of cpus to utilize)",
	action="store",
	dest="cpu")

parser.add_argument("-gs", "-guide_start",
	required=False,
	default=0,
	help="position of guide start in read2",
	action="store",
	dest="guide_start")

parser.add_argument("-bl", "-barcode_length",
	required=False,
	default=31,
	help="length of barcode (default=31)", 
	action="store", 
	dest="barcode_length")

parser.add_argument("-rc", "-read_cutoff",
	required=False,
	default=0,
	help="read count cutoff for barcodes to keep (default=0)",
	action="store",
	dest="read_cutoff")

parser.add_argument("-cc", "-cluster_cutoff",
	required=False,
	default=10,
	help="clusteringg count cutoff",
	action="store",
	dest="cluster_cutoff")

parser.add_argument("-bq", "-bquality",
	required=False,
	default=53,
	help="ascii quality score cutoff for barcode (default=53)",
	action="store",
	dest="barcode_quality")

parser.add_argument("-gdq", "-gdquality",
	required=False,
	default=55,
	help="ascii quality score cutoff for guide-donor (default=55)",
	action="store",
	dest="guide_donor_quality")

parser.add_argument("-psl", "-primingseqlength",
	required=False,
	default=20,
	help="length of priming sequence in bps to exclude from analysis (default=20)",
	action="store",
	dest="priming_seq_length")

parser.add_argument("-bcd", "-barcode_cluster_distance",
	required=False,
	default=4,
	help="maximum hamming distance between barcodes to permit mergeing (default=4)",
	action="store",
	dest="barcode_cluster_distance")

parser.add_argument("-gcd", "-guide_cluster_distance",
	required=False,
	default=2,
	help="maximum hamming distance between guides to permit mergeing (default=2)",
	action="store",
	dest="guide_cluster_distance")

parser.add_argument("-dcd", "-donor_cluster_distance",
	required=False,
	default=8,
	help="maximum hamming distance between donors to permit mergeing (default=8)",
	action="store",
	dest="donor_cluster_distance")


BCD = int(args.barcode_cluster_distance)
GCD = int(args.guide_cluster_distance)
DCD = int(args.donor_cluster_distance)

args = parser.parse_args()

# Constant variables
PATH = os.path.realpath(args.path)
TMP = PATH + "/.tmp/"

BARCODE_LENGTH = int(args.barcode_length) # size of barcode
READ_COUNT_CUTOFF = int(args.read_cutoff) # cutoff for reads to exclude, default is zero
BARCODE_QUALITY_CUTOFF = int(args.barcode_quality) # cutoff for read quality
GUIDE_DONOR_QUALITY_CUTOFF = int(args.guide_donor_quality) # cutoff for read quality
TOTAL_SEGMENTS = int(args.cpu) # total number of segments to run processes at (multithreading)
DESIGN_FILE = args.design_file
DATABASE = args.database
LIBRARY_NAME = PATH.split("/")[-1]
CLUSTER_CUTOFF = args.cluster_cutoff
GUIDE_START = int(args.guide_start)
R2_FILTER = args.r2_filter
PRIMING_SEQ_LENGTH = args.priming_seq_length

def load_reads(files):
	lines = []
	for file in files:
		lines.extend(gzip.open(file).readlines())

	sequence = [lines[r] for r in range(1, len(lines), 4)] # read sequences
	sequence = [l.decode('utf-8').replace("\n","") for l in sequence]

	quality = [lines[r] for r in range(3, len(lines), 4)] # read quality scores
	quality = [l.decode('utf-8').replace("\n","") for l in quality]

	return sequence, quality

# Assert correct structure of the library design directory
print(str(datetime.datetime.now()) + " Checking for required precursor files.")
try:
	assert os.path.isdir(PATH)
	assert os.path.isdir(PATH + "/sequencing/read1")
	assert os.path.isdir(PATH + "/sequencing/read2")

except AssertionError:
    print(str(datetime.datetime.now()) + " Incorrect layout for the library design directory.")
    raise

if not (os.path.isdir(TMP)):
	os.mkdir(TMP)

print(str(datetime.datetime.now()) + " Checking to see if any intermediate files already exist.")
temp_files = os.listdir(TMP)

pickle_files = ["read_count_dict"] + [str(seg) + "-" + str(TOTAL_SEGMENTS) + ".R1_dict" for seg in range(0, TOTAL_SEGMENTS)] + [str(seg) + "-" + str(TOTAL_SEGMENTS) + ".R2_dict" for seg in range(0, TOTAL_SEGMENTS)]

skip_pickling = True
for i in pickle_files:
	if i not in temp_files:
		skip_pickling = False
		break;

count_dict_pickle = TMP + "read_count_dict"

if skip_pickling:
	print(str(datetime.datetime.now()) + " Pickled barcode to read mappings already present. Skipping!")
else:
	# Collect sequencing files.
	print(str(datetime.datetime.now()) + " Creating temporary pickle files for reads sorted by barcode, for " + str(TOTAL_SEGMENTS) + " intended threads.")
	read1_files = [i for i in os.listdir(PATH + "/sequencing/read1") if ("fastq.gz" in i)]
	read2_files = [i for i in os.listdir(PATH + "/sequencing/read2") if ("fastq.gz" in i)]

	# Load reads from specified files.
	print(str(datetime.datetime.now()) + " Reading read1 files into memory.")
	read1_sequence, read1_quality = load_reads([PATH + "/sequencing/read1/" + f for f in read1_files])

	print(str(datetime.datetime.now()) + " Reading read2 files into memory.")
	read2_sequence, read2_quality = load_reads([PATH + "/sequencing/read2/" + r for r in read2_files])

	if (R2_FILTER != None): # Optionally filter out reads lacking a common priming sequence.
		inital_len = len(read1_sequence)
		print(str(datetime.datetime.now()) + " Remove reads missing common priming sequencing in R2.")
		read1_sequence, read1_quality, read2_sequence, read2_quality = zip(*[(a, b, c, d) for a, b, c, d
			in zip(read1_sequence, read1_quality, read2_sequence, read2_quality) if c[:GUIDE_START] == R2_FILTER])

		print(str(datetime.datetime.now()) + " " + str(inital_len-len(read1_sequence)) + " reads removed.")
		
	if (GUIDE_START != 0): # Optionally trim reads.
		print(str(datetime.datetime.now()) + " Trimming read2 files.")
		read2_sequence = [i[GUIDE_START:] for i in read2_sequence]
		read2_quality = [i[GUIDE_START:] for i in read2_quality]

	print(str(datetime.datetime.now()) + " Calculating quality scores for read1 files.")
	barcode_quality_scores = []
	read1_guide_donor_quality_scores = []

	counter = 0
	for line in read1_quality: # Get ASCII quality score for read1 and for the barcode portion alone
	        barcode_scores = [ord(i) for i in line[:BARCODE_LENGTH]]
	        barcode_quality_scores.append(np.mean(barcode_scores))

	        gd_scores = [ord(i) for i in line[BARCODE_LENGTH:]]
	        read1_guide_donor_quality_scores.append(np.mean(gd_scores))

	        counter +=1
	        if(counter % 1000000 == 0):
	        	print(str(datetime.datetime.now()) + " " + str(counter) + " reads and counting.")

	print(str(datetime.datetime.now()) + " Calculating quality score for read2 files.")
	read2_guide_donor_quality_scores = []

	counter = 0
	for line in read2_quality: # Get ASCII quality score for read2
	        scores = [ord(i) for i in line]
	        read2_guide_donor_quality_scores.append(np.mean(scores))

	        counter +=1
	        if(counter % 1000000 == 0):
	        	print(str(datetime.datetime.now()) + " " + str(counter) + " reads and counting.")

	print(str(datetime.datetime.now()) + " Filtering out low-quality reads.")
	initial_number_reads = len(read1_sequence)
	read1_sequence, read2_sequence, barcodes = zip(*[(f[BARCODE_LENGTH+PRIMING_SEQ_LENGTH:], r, f[:BARCODE_LENGTH]) for f, r, fscore, fscore2, rscore
	        in zip(read1_sequence, read2_sequence, barcode_quality_scores, read1_guide_donor_quality_scores, read2_guide_donor_quality_scores) 
	        if (fscore >= BARCODE_QUALITY_CUTOFF) and (fscore2 >= GUIDE_DONOR_QUALITY_CUTOFF) and (rscore >= GUIDE_DONOR_QUALITY_CUTOFF)]) #f[:BARCODE_LENGTH-20]

	print(str(datetime.datetime.now()) + " " + str(initial_number_reads - len(read1_sequence)) + " of " + str(initial_number_reads) + " reads excluded.")

	if (READ_COUNT_CUTOFF != 0): # Optional step to remove low read barcodes from annotations.
		print(str(datetime.datetime.now()) + " Read count cutoff option selected, removing low read barcodes from analysis.")
		barcodes_to_keep = [key for key, count in Counter(barcodes).items() if count >= READ_COUNT_CUTOFF]
		keep_dict = {g: True for g in barcodes_to_keep}
		read1_sequence, read2_sequence, barcodes = zip(*[(f, r, b) for f, r, b 
			in zip(read1_sequence, read2_sequence, barcodes) if b in keep_dict])

	# Store barcode read count dictionary for later use. 
	count_dict = dict(Counter(barcodes))
	print(str(datetime.datetime.now()) + " " + str(len(count_dict)) + " barcodes in total.")

	pickle_out = open(count_dict_pickle, "wb")
	pickle.dump(count_dict, pickle_out, protocol=2)
	pickle_out.close()

	# Divide up barcodes into specified number of segments for parallel analysis.
	LENGTH = len(set(barcodes))
	total_segments = int(TOTAL_SEGMENTS)

	barcode_list = list(set(barcodes))
	for segment in range(0, total_segments):
		start = int((LENGTH/total_segments)*segment) # determine start and end position of segment.
		if (segment+1 == total_segments):
			sub_barcodes_set = barcode_list[start:]
		else:
			stop = int((LENGTH/total_segments)*(segment+1))
			sub_barcodes_set = barcode_list[start:stop]
		sub_barcodes_dict = {b: True for b in sub_barcodes_set}

		sub_read1, sub_read2, sub_barcodes = zip(*[(f, r, b) for f, r, b 
			in zip(read1_sequence, read2_sequence, barcodes) if b in sub_barcodes_dict])

		R1_dict, R2_dict = {}, {} # store reads by barcode into R1 and R2 dictionaries.
		for f, r, b in zip(sub_read1, sub_read2, sub_barcodes):
			if (b not in R1_dict) and (b not in R2_dict):
				R1_dict[b] = [f]
				R2_dict[b] = [r]
			else:
				R1_dict[b].append(f)
				R2_dict[b].append(r)

		seg_unit = str(segment) + "-" + str(total_segments)

		pickle_out = open(TMP + seg_unit + ".R1_dict", "wb")
		pickle.dump(R1_dict, pickle_out, protocol=2)
		pickle_out.close()

		pickle_out = open(TMP + seg_unit + ".R2_dict", "wb")
		pickle.dump(R2_dict, pickle_out, protocol=2)
		pickle_out.close()


consensus_files = [str(seg) + "-" + str(TOTAL_SEGMENTS) + '_consensus_queries.fasta' for seg in range(0, TOTAL_SEGMENTS)]
skip_consensus_finding = True
for i in consensus_files:
	if i not in temp_files:
		skip_consensus_finding = False
		break;

if skip_consensus_finding:
	print(str(datetime.datetime.now()) + " Consensus sequences for guide-donor pairs already determined. Skipping!")
else:
	print(str(datetime.datetime.now()) + " Beginning threaded executions of consensus finding process.")
	commands = []
	for i in range(0, TOTAL_SEGMENTS):
		seg_unit = str(i) + "-" + str(TOTAL_SEGMENTS)
		f1_name = TMP + seg_unit + ".R1_dict"
		f2_name = TMP + seg_unit + ".R2_dict"

		commands.append(subprocess.Popen("python findConsensus.py -r1 " + f1_name + 
			" -r2 " + f2_name + 
			" -o " + TMP + 
			" -s " + str(i) + 
			" -ts " + str(TOTAL_SEGMENTS), shell=True))

	exit_codes = [p.wait() for p in commands]
	print(str(datetime.datetime.now()) + " Consensus finding complete.")

blastn_files = [str(seg) + "-" + str(TOTAL_SEGMENTS) + "_blast.csv" for seg in range(0, TOTAL_SEGMENTS)]
skip_blastn = True
for i in blastn_files:
	if i not in temp_files:
		skip_blastn = False
		break;

if skip_blastn:
	print(str(datetime.datetime.now()) + " Blastn alignments already completed. Skipping!")
else:
	print(str(datetime.datetime.now()) + " Beginning threaded executions of blastn.")
	commands = []

	for j in range(0, TOTAL_SEGMENTS):
		consensus_file = TMP + str(j) + "-" + str(TOTAL_SEGMENTS) + '_consensus_queries.fasta'
		commands.append(subprocess.Popen('blastn -query ' + consensus_file + 
			' -db ' + DATABASE + 
			' -outfmt "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" -out ' + 
			consensus_file.split("consensus")[0] + 'blast.csv -num_alignments 1', shell=True))

	exit_codes = [p.wait() for p in commands]

blastn_files = [str(seg) + "-" + str(TOTAL_SEGMENTS) + "_blast_mapped.csv" for seg in range(0, TOTAL_SEGMENTS)]
skip_blastn_alignment = True
for i in blastn_files:
	if i not in temp_files:
		skip_blastn_alignment = False
		break;

if skip_blastn_alignment:
	print(str(datetime.datetime.now()) + " Blastn alignments annotations already completed. Skipping!")
else:
	print(str(datetime.datetime.now()) + " Beginning threaded executions of blastn alignment annotations to design file.")
	commands = []

	for j in range(0, TOTAL_SEGMENTS):
		consensus_file = TMP + str(j) + "-" + str(TOTAL_SEGMENTS) + '_consensus_queries.fasta'
		commands.append(subprocess.Popen("python runBlastn.py -in " + consensus_file + 
			" -db " + DATABASE + 
			" -p " + DESIGN_FILE, shell=True))

	exit_codes = [p.wait() for p in commands]

if "unclustered_reference_table.csv" in temp_files:
	print(str(datetime.datetime.now()) + " Guide-donor statuses already annotated. Skipping!")
else:
	blast_files = []
	for k in range(0, TOTAL_SEGMENTS):
		blast_files.append(TMP + str(k) + "-" + str(TOTAL_SEGMENTS) + "_blast_mapped.csv")

	frames = []
	for i in range(0, TOTAL_SEGMENTS):
		frames.append(pd.read_csv(blast_files[i]))
	result = pd.concat(frames)

	pickle_in = open(count_dict_pickle, "rb")
	count_dict = pickle.load(pickle_in)
	count_df = pd.DataFrame.from_dict(count_dict, orient='index')
	count_df.columns = ['read_count']

	result['qseqid'] = [i.split("_")[0] for i in list(result['qseqid'])]
	result = pd.merge(result, count_df, left_on="qseqid", right_index=True, how="left")
	result['guide'] = [i[:20] for i in result['qseq']]
	result['donor'] = [i[BARCODE_LENGTH:] for i in result['qseq']]

	inital_len = len(result)
	#result = result[(result['length'] >= 125) & (result['length'] <= 155)]

	print(str(datetime.datetime.now()) + " Annotating guide and donor status.")
	#print(str(datetime.datetime.now()) + " Excluding " + str(inital_len - len(result)) + " of " + str(inital_len) + " guide-donor pairs with abnormal genomic mapping lengths.")

	# Annotate guide and donor status.
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
	result['library'] = LIBRARY_NAME

	print(str(datetime.datetime.now()) + " Guide status: " + str(dict(result['guide_status'].value_counts())))
	print(str(datetime.datetime.now()) + " Donor status: " + str(dict(result['donor_status'].value_counts())))
	result.to_csv(TMP + "unclustered_reference_table.csv", index=False)

if "unclustered_reference_table.csv" not in temp_files:	
	print(str(datetime.datetime.now()) + " Clustering barcodes in parallel.")
	commands = []

	for i in range(0, TOTAL_SEGMENTS):
		consensus_file = TMP + str(i) + "-" + str(TOTAL_SEGMENTS) + '_consensus_queries.fasta'
		commands.append(subprocess.Popen("python combineRefAlign.py -in " + TMP + "unclustered_reference_table.csv " + 
			" -cutoff " + str(CLUSTER_CUTOFF) + 
			" -segment " + str(i) + 
			" -total_segments " + str(TOTAL_SEGMENTS) +  
			"-bcd " + str(BCD) +
			"-gcd " + str(GCD) + 
			"-dcd " + str(DCD),shell=True))

	exit_codes = [p.wait() for p in commands]

files = [TMP + str(segment) + ".cluster_dict" for segment in range(0, TOTAL_SEGMENTS)]
ff = TMP + "unclustered_reference_table.csv"

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
#ref_table = ref_table[~ref_table.index.duplicated(keep='first')]
pre_clustering_length = len(ref_table)
ref_table['read_count'] = ref_table.groupby(['cluster_id'])['read_count'].transform('sum')
ref_table.drop_duplicates('cluster_id', inplace=True)
ref_table.to_csv(PATH + "/" + PATH.split("/")[-1] + "_reference_table.csv")

print(str(datetime.datetime.now()) + " Clustering from " + str(pre_clustering_length) + " to " + str(len(ref_table)) + " barcodes.")
print(str(datetime.datetime.now()) + " Done!")

