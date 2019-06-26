# MAGESTIC library mapping and annotation pipeline

Aim: to generically process pair-end reads of barcode-guide-donor cassettes and create a library reference table mapping barcodes to features.

### Introduction

These scripts are currently designed to use MAGESTIC library paired-end sequencing data for:
1. mapping barcodes to designed guide-donor pairs,
2. identifying synthesis mutation errors in the guide-donor sequence, and
3. cluster barcodes using sequence similarity and relative abundance criteria.

The library design used as a model during pipeline development was:
- 25 nt pseudo-randomized barcode
- 110 nt donor sequence
- 11 nt BspQI restriction enzyme cut site
- 20 nt guide sequence

The paired-end design used as a model during pipeline development is as follows:
- 151 bp "read1" with the barcode + primer sequence + donor sequence
- 151 bp "read2" with the guide + BspQI cut site + partial donor sequence

Although these were the intended designed formats, the pipeline is somewhat adaptable for future library designs using optional arguments.

### Basic Command

Given that the necessary dependencies are present, the pipeline can be run by executing run.py from the command line as follows.

```
usage: run.py [-h] -d PATH -r DESIGN_FILE -db DATABASE [-c CPU]
              [-gs GUIDE_START] [-bl BARCODE_LENGTH] [-r2f R2_FILTER]
              [-rc READ_CUTOFF] [-cc CLUSTER_CUTOFF] [-bq BARCODE_QUALITY]
              [-gdq GUIDE_DONOR_QUALITY]

optional arguments:
  -h, --help            show this help message and exit
  -d PATH, -dir PATH    name of library design directory
  -r DESIGN_FILE, -ref DESIGN_FILE
                        reference guide-donor design file
  -db DATABASE, -database DATABASE
                        blastn database to use for mapping
  -c CPU, -cpu CPU      the number of partitions to split data into (i.e. num
                        of cpus to utilize)
  -gs GUIDE_START, -guide_start GUIDE_START
                        position of guide start in read2
  -bl BARCODE_LENGTH, -barcode_length BARCODE_LENGTH
                        length of barcode (default=31)
  -r2f R2_FILTER, -read2_filter R2_FILTER
  -rc READ_CUTOFF, -read_cutoff READ_CUTOFF
                        read count cutoff for barcodes to keep (default=0)
  -cc CLUSTER_CUTOFF, -cluster_cutoff CLUSTER_CUTOFF
                        clusteringg count cutoff
  -bq BARCODE_QUALITY, -bquality BARCODE_QUALITY
                        ascii quality score cutoff for barcode (default=53)
  -gdq GUIDE_DONOR_QUALITY, -gdquality GUIDE_DONOR_QUALITY
                        ascii quality score cutoff for guide-donor
                        (default=55)
```

For example, for the REDI validation pipelines, the following command could be used:
```
python run.py -d V313_sample1_step1_20190605 
	-r fasta_files/removed_subpool_priming_sequence/RM11_PAM_disruption1_subpool_1.fasta.gz 
	-db blastn_databases/remove_subpool_priming_sequence/RM11_PAM_disruption1_subpool_1_mod 
	-c 6 -gs 28 -bl 25 > sample1.out &
```

### Dependencies

In order to execute this pipeline, an installed copy of BLASTN is required and must be used to create a database mapping designed features to their corresponding guide-donor sequences. This can be done by following the steps below:

1. Create a library design FASTA file with the following format:
	```
	>chrX_18819_TCA_to_GCC_PAM_18817_linked_2.1:PAM_disruption_or_1_to_5_bp_from_PAM_disruption
	GCTATAACGGGCCAAGTGACgtttgaagagcATGCAATAATGGGGTCTGCTTGCAGATGGTATAATCTGCTATAACGGGCCAAGGGCCAGGCACGCTGTTGCAGAACATACCATCTTATATCATGTAGGTTGTGCAAACA
	>chrX_18819_TCA_to_GCC_PAM_18823_+_linked_3.1:PAM_disruption_or_1_to_5_bp_from_PAM_disruption
	GCAACAGCGTGCCTGTCACTgtttgaagagcTGTTTGCACAACCTACATGATATAAGATGGTATGTTCTGCAACAGCGTGCCTGGCCCTTGGCCCGTTATAGCAGATTATACCATCTGCAAGCAGACCCCATTATTGCAT
	>chrX_24374_CTGC_to_TTGT_PAM_24378_+_linked_2.1:PAM_disruption_or_1_to_5_bp_from_PAM_disruption
	TCTTCCTTCCCATGGCTGCAgtttgaagagcTTTATGTTAAATGTGCTCGTTACTTCAACAAGTAAAGTCTTCCTTCCCATGGTTGTAAGGTGCGGGAATCACCATTTGAGTTTGCATAGTTGCCAAAAAACATGGTGTA
	```
where the sequence headers are the feature descriptors and the sequence itself is the guide-donor sequence expected from the paired-end reads.

2. A library database must the be generated using the library design file and the ```makeblastdb``` command from BLASTN.
	```
	./makeblastdb -in [LIBRARY_DESIGN_FILE] 
			-dbtype nucl 
			-out [DATABASE_NAME]
	```

3. Finally, either a direct or symbolic link to the ```runblastdb``` package must be in the working directory, so that BLASTN alignments can be performed using your generated library database and the alignment sequences.

For each new library developed, this process will need to be repeated, and the Snakefile described below will need to be updated to include information relevant to that specific library.

### Sequencing Data Format

Sequencing data must be in the gzipped fastq format according to the following directory structure.

path/V313_sample1_step1_20190605
path/V313_sample1_step1_20190605/sequencing
path/V313_sample1_step1_20190605/sequencing/read1
path/V313_sample1_step1_20190605/sequencing/read2

This means that for a command of run.py where -d V313_sample1_step1_20190605, the sequencing data should be allocated into appropriate read1 and read2 files. 

During script execution, a temporary directory will be created to store intermediate files.
path/V313_sample1_step1_20190605/.tmp/

### Scripts

Below is an overview of the pipeline, starting with an explanation of the arguments to be provided upon execution of run.py