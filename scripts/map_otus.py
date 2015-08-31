#! /usr/bin/env python3

# Run usearch to map sequences to OTUs

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-06-05

#  **** NOTE ****
#  This version assumes the unique sequences are still in FASTA files in the working
#  directory of the remove_duplicates script.  The default directory name is 'uniq'
#  but an alternative can be specified with --directory.
#
#  TBD: don't have a default directory; if no directory is specified write the sequences
#  in the uniq table to the working directory for this script.
#  ***************

import sqlite3
import argparse
import os
import os.path
import sys

from common import *
from FASTQ import *

# Define filenames and patterns (may be referenced by more than one function)

result_file_pattern = 'readmap.{}.uc'
input_file_pattern = 'unique.{}.fasta'
ref_db_file = 'otus.fasta'

###
# Make the reference "database" (FASTA file) from non-chimeric OTU seeds

fetch_sequences = "SELECT cluster_id, sequence FROM (SELECT cluster_id, count(members.name) AS n, sequence FROM clusters LEFT JOIN members USING (cluster_id) GROUP BY cluster_id) JOIN chimeras USING (cluster_id) WHERE n > 0 AND chimeric = 'N'"

def make_reference_db(db, args):
    ff = open(os.path.join(args.workspace, ref_db_file), 'w')
    record_metadata(db, 'query', fetch_sequences)
    for otu_id, sequence in db.execute(fetch_sequences).fetchall():
        print('>OTU_{}'.format(otu_id), file=ff)
        print(sequence, file=ff)
    ff.close()

###
# Run the app

def run_usearch_global(sid, args):
    cmnd = 'usearch -usearch_global '
    cmnd += os.path.join(args.directory, input_file_pattern.format(sid))
    cmnd += ' -db ' + os.path.join(args.workspace, ref_db_file)
    cmnd += ' -strand plus'
    cmnd += ' -id 0.97'
    cmnd += ' -uc ' + os.path.join(args.workspace, result_file_pattern.format(sid))
    print(cmnd)
    record_metadata(db, 'exec', cmnd, commit=True)
    res = os.system(cmnd)

# Populate the table 

insert_record = 'INSERT INTO otus (otu_id, sample_id, count) VALUES (?,?,?)'
fetch_counts = 'SELECT defline, n FROM panda JOIN uniq USING (panda_id) where uniq.sample_id = {}'

def import_results(db, args, sid):
    cmap = defline_map(db, sid)            # map deflines to number of times seq found in this sample
    count = { }
    for line in open(os.path.join(args.workspace, result_file_pattern.format(sid))):
        res = line.split('\t')
        otu = res[-1].strip()
        otu_id = 0 if otu == '*' else int(otu.split('_')[-1])
        defline = FASTQ.parse_defline(res[8])
        count.setdefault(otu_id, 0)
        count[otu_id] += cmap[defline]
    for otu_id in sorted(count.keys()):
        db.execute(insert_record, (otu_id, sid, count[otu_id]))

def defline_map(db, sid):
    dm = { }
    query = fetch_counts.format(sid)
    for defline, n in db.execute(query):
        dm[defline] = n
    return dm
    
###
# Top level function: initialize the workspace directory, run the app

def map_otus(db, args):
    init_workspace(args)
    make_reference_db(db,args)
    for row in sample_list(db, args):
        sid = row[0]
        run_usearch_global(sid, args)
        import_results(db, args, sid)

###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    
    args = init_api(
        desc = "Run usearch to map merged sequences to one of the inferred clusters.",
        specs = [
            ('workspace',    { 'metavar': 'dir', 'help' : 'working directory', 'default' : 'map' } ),
            ('directory',    { 'metavar': 'dir', 'help' : 'name of directory containing unique sequences', 'default' : 'uniq' } ),
            ('sample',       { 'metavar': 'id', 'help' : 'process sequences from this sample only'} )
        ]
    )
        
    db = sqlite3.connect(args.dbname)
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))
    
    try:
        otu_spec = [('sample_id', 'foreign', 'samples'), ('count', 'INTEGER')]
        init_table(db, 'otus', 'otu_id', otu_spec, args.force, has_primary=False)
    except Exception as err:
        print('Error while initializing output tables:', err)
        argparse.ArgumentParser.exit(1, 'Script aborted')

    map_otus(db, args)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )
    record_metadata(db, 'end', '')

    db.commit()
