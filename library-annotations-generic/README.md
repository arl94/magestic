
# MAGESTIC library mapping and annotation pipeline

This pipeline is intended to map pair-end reads of barcode-guide-donor cassettes to the designed variants and create a library reference table with annotated information on the barcode-guide-donor pairs.

### Introduction

These scripts are currently designed to use MAGESTIC library paired-end sequencing data for:
1. mapping barcodes to designed guide-donor pairs,
2. identifying synthesis mutation errors in the guide-donor sequence, and
3. clustering barcodes using sequence similarity and relative abundance criteria.

The pipeline accepts a variety of guide/donor/barcode sizes and positions (i.e., the format of the library can vary). However, there are some limitations for how customizable the mapping and annotation process with the current command arguments.

Feel free to adjust the parameters in the code to satisfy your needs.

### Dependencies

In order to execute this pipeline, an installed copy of BLASTN is required and must be used to create a database mapping designed features to their corresponding guide-donor sequences. The necessary package for BLASTN is included with this GitHub repository.

Creating the necessary BLASTN database is done by following the steps below:

1. Produce a library design FASTA file with the following format:
	```
	>chrX_18819_TCA_to_GCC_PAM_18817_linked_2.1:PAM_disruption_or_1_to_5_bp_from_PAM_disruption
	GCTATAACGGGCCAAGTGACgtttgaagagcATGCAATAATGGGGTCTGCTTGCAGATGGTATAATCTGCTATAACGGGCCAAGGGCCAGGCACGCTGTTGCAGAACATACCATCTTATATCATGTAGGTTGTGCAAACA
	>chrX_18819_TCA_to_GCC_PAM_18823_+_linked_3.1:PAM_disruption_or_1_to_5_bp_from_PAM_disruption
	GCAACAGCGTGCCTGTCACTgtttgaagagcTGTTTGCACAACCTACATGATATAAGATGGTATGTTCTGCAACAGCGTGCCTGGCCCTTGGCCCGTTATAGCAGATTATACCATCTGCAAGCAGACCCCATTATTGCAT
	>chrX_24374_CTGC_to_TTGT_PAM_24378_+_linked_2.1:PAM_disruption_or_1_to_5_bp_from_PAM_disruption
	TCTTCCTTCCCATGGCTGCAgtttgaagagcTTTATGTTAAATGTGCTCGTTACTTCAACAAGTAAAGTCTTCCTTCCCATGGTTGTAAGGTGCGGGAATCACCATTTGAGTTTGCATAGTTGCCAAAAAACATGGTGTA
	```
where the sequence headers are the feature descriptors and the sequence itself is the guide-donor sequence expected from the paired-end reads. With MAGESTIC, this is the library design file specifying which guide-donor pairs have been introduced into the library.

2. Generate a library database using the library design file from the prior step and the ```makeblastdb``` command from BLASTN.
	```
	./makeblastdb -in [LIBRARY_DESIGN_FILE] 
			-dbtype nucl 
			-out [DATABASE_NAME]
	```

3. Finally, either a direct or symbolic link to the ```runblastdb``` package must be in the working directory, so that BLASTN alignments can be performed using your generated library database and the alignment sequences. The script for ```runblastdb``` is included with this repository.

For each new library developed, a new library database must be generated. This is the file which aligned guide-donor pairs are mapped to from library sequencing information, and thus it should be as accurate as possible.

### Basic Command

Given that the necessary dependencies are present, the pipeline can be run by executing run.py from the command line as follows.

```
usage: run.py [-h] -d PATH -r DESIGN_FILE -db DATABASE [-c CPU]
              [-gs GUIDE_START] [-bl BARCODE_LENGTH] [-r2f R2_FILTER]
              [-rc READ_CUTOFF] [-cc CLUSTER_CUTOFF] [-bq BARCODE_QUALITY]
              [-gdq GUIDE_DONOR_QUALITY]

mandatory arguments:
-d PATH, -dir PATH
	name of library design directory
-r DESIGN_FILE, -ref DESIGN_FILE
	reference guide-donor design file
-db DATABASE, -database DATABASE
	blastn database to use for mapping

optional arguments:
-c CPU, -cpu CPU
	the number of partitions to split data into (i.e. num of cpus to utilize)
-gs GUIDE_START, -guide_start GUIDE_START
	position of guide start in read2
-bl BARCODE_LENGTH, -barcode_length BARCODE_LENGTH
	length of barcode (default=31)
-rc READ_CUTOFF, -read_cutoff READ_CUTOFF
	read count cutoff for barcodes to keep (default=0)
-cc CLUSTER_CUTOFF, -cluster_cutoff CLUSTER_CUTOFF
	clustering count cutoff (only barcodes below this threshold will be collapsed, and only into barcodes above this threshold) (default=10)
-bq BARCODE_QUALITY, -bquality BARCODE_QUALITY
	ascii quality score cutoff for barcode (default=53)
-gdq GUIDE_DONOR_QUALITY, -gdquality GUIDE_DONOR_QUALITY
	ascii quality score cutoff for guide-donor (default=55)
```

For example, for the REDI validation pipelines, the following command could be used:
```
python run.py -d V313_sample1_step1_20190605 
	-r fasta_files/removed_subpool_priming_sequence/RM11_PAM_disruption1_subpool_1.fasta.gz 
	-db blastn_databases/remove_subpool_priming_sequence/RM11_PAM_disruption1_subpool_1_mod 
	-c 6 -gs 28 -bl 25 > sample1.out &
```

In this above example, the BLASTN database that has been created is called "RM11_PAM_disruption1_subpool1_1_mod", the guide sequence starts at nucleotide position 28 in read 2, and the barcode length is 25 nts.

Note that the guide must be in read 2. The general format of the reads should be:
Read 1: barcode - donor
Read 2: priming sequence (optional) - guide - BspQI cut site - donor

### Sequencing Data Format

Sequencing data must be in the fastq.gz (gzip) format according to the following directory structure:

	path/V313_sample1_step1_20190605
	path/V313_sample1_step1_20190605/sequencing
	path/V313_sample1_step1_20190605/sequencing/read1
	path/V313_sample1_step1_20190605/sequencing/read2

This means that for a command of run.py where -d V313_sample1_step1_20190605, the sequencing data should be allocated into appropriate read1 and read2 files. 

During script execution, a temporary directory will be created to store intermediate files.
	
	path/V313_sample1_step1_20190605/.tmp/

The ```.tmp/``` folder stores intermediate files. These can be deleted upon successful completion of the pipeline. If the pipeline fails or prematurely terminates, and the command is restarted, it will pick up midway in the pipeline using these temporary files.

### Overview

Below is an overview of the pipeline.

1) Quality check: read1 and read2 files are filtered using the quality scores from sequencing the cutoffs indicated with the ```-bq``` and ```-gdq``` arguments.
2) Barcode pickling: read1 and read2 are pickled into dictionaries mapping a barcode to its corresponding reads.
3) Partitioning: the barcodes are partitioned using the ```-cpu``` argument so that later steps can be completed in a multi-threaded manner.
4) Consensus finding: for each barcode, a read1 and read2 consensus is found individually; then, the read1 and read2 consensus are aligned to form a single, complete consensus.
5) Feature mapping: the consensus sequence for each barcode is mapped to its most similar guide-donor feature from the indicated BLASTN database (see: Dependencies).
6) Feature annotation: the barcode-feature mappings are then annotated with guide_status, donor_status, and more (note: these annotations will quickly be unreliable if the format of read1 and/or read2 change significantly.
7) Clustering: lowly abundant barcodes are collapsed into barcodes of greater abundance, based on the cutoffs specified in ```-cc```, ```-bcd```, ```-gcd```, and ```-dcd``` parameters.
8) The final reference table is outputted in the directory and program terminates.