META 16S rRNA Analysis Pipeline

John Conery
META Center for Systems Biology
University of Oregon

http://meta.uoregon.edu/

This project contains a set of Python scripts for managing applications 
used in an analysis pipeline for 16S rRNA genes.  The software was 
originally developed for microbial communities sampled from zebrafish 
but should be easily adapted for other microbial communities.

The distinguishing feature of this pipeline is the use of SQLite to 
manage the workflow.  After the first few stages reduce the data via
quality filtering and dereplication the sequence data and all the
derived information is saved in a single SQLite database file.  

Pipeline Stages
---------------
db_setup.py             initialize a new project
import_reads.py         specify locations of .fastq files
assemble_pairs.py       combine paired ends into single sequences
remove_duplicates.py    dereplication
map_reference.py        find known sequences to seed clusters
form_otus.py            clustering
filter_chimeras.py      identify chimeric sequences
map_otus.py             assign sequences to clusters
classify.py             taxonomic classification of clusters

Data Analysis
-------------
abundance.py            group data by taxonomic classifications
print_abundance.py      print abundances in CSV format
summarize.py            print table sizes

Scripts for Comparing Databases
-------------------------------
compare_clusters.py     find clusters identified in two databases
compare_members.py      find sequences assigned to different clusters
compare_otus.py         compare taxonomic classes in two databases

Python Modules and Other Scripts
--------------------------------
FASTA.py                class definition for FASTA sequences
FASTQ.py                class definition for FASTQ sequences
fstats.py               print statistics about FASTA or FASTQ file
common.py               functions used by most scripts
config.py               paths to external applications
gref.py                 print sequences matching a pattern
make_script.py          generate a script to run the pipeline
print_as_fastq.py       print a sequence table in FASTQ format

Deprecated
----------
quality_filter.py       subsumed by assemble_pairs

