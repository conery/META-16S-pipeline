#! /usr/bin/env python3

# Run usearch to filter chimeras using a reference database

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-06-05

import sqlite3
import argparse
import os
import os.path
import re
import sys

from common import *

# Define filenames (may be referenced by more than one function)

result_file = 'results.uc'
input_file = 'clusters.fasta'

###
# Print the sequences to check; they're the centroids in the clusters saved in the clusters table


def print_sequences(db, args):
    fetch_clusters = 'SELECT name, sequence FROM clusters'
    record_metadata(db, 'query', fetch_clusters)
    ff = open(os.path.join(args.workspace, input_file), 'w')
    for name, sequence in db.execute(fetch_clusters):
        print('>{}'.format(name), file=ff)
        print(sequence, file=ff)    
    ff.close()

###
# Run the app

def run_uchime_ref(args):
    cmnd = 'usearch -uchime_ref '
    cmnd += os.path.join(args.workspace, input_file)
    cmnd += ' -db ' + args.reference
    cmnd += ' -strand plus'
    cmnd += ' -uchimeout ' + os.path.join(args.workspace, result_file)
    print(cmnd)
    record_metadata(db, 'exec', cmnd, commit=True)
    res = os.system(cmnd)
    
###
# Parse the output file, save references to chimeric sequences in a new table.

insert_record = 'INSERT INTO chimeras (cluster_id, chimeric) VALUES (?,?)'

def import_results(db, args):
    cmap = make_cmap(db)
    for line in open(os.path.join(args.workspace, result_file)):
        res = line.split('\t')
        defspec = res[1]
        chimeric = res[-1].strip()
        db.execute(insert_record, (cmap[defspec], chimeric))

def make_cmap(db):
    cmap = { }
    for cid, name in db.execute('SELECT cluster_id, name FROM clusters'):
        cmap[name] = cid
    return cmap

###
# Top level function: initialize the workspace directory, run the app

def filter_chimeras(db, args):
    init_workspace(args)
    print_sequences(db, args)
    run_uchime_ref(args)
    import_results(db, args)

###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    
    args = init_api(
        desc = "Run usearch to filter chimeras using a reference sequence.",
        specs = [
            ('workspace',    { 'metavar': 'dir', 'help' : 'working directory', 'default' : 'chimeras' } ),
            ('reference',    { 'metavar': 'fn', 'help' : 'FASTA file containing reference sequences', 'default' : path_to_resource('gold.fa') } ),
        ]
    )
        
    db = sqlite3.connect(args.dbname)
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))

    try:
        chimera_spec = [('cluster_id', 'foreign', 'clusters'), ('chimeric', 'CHAR(1)')]
        init_table(db, 'chimeras', 'chimera_id', chimera_spec, args.force)
    except Exception as err:
        print('Error while initializing output tables:', err)
        argparse.ArgumentParser.exit(1, 'Script aborted')

    filter_chimeras(db, args)
    record_metadata(db, 'end', '')

    db.commit()
