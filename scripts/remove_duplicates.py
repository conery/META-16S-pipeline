#! /usr/bin/env python3

# Run cd-hit-dup to find unique assembled reads.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-05-07

#  **** NOTE ****
#  This version assumes assemble_pairs was run with the --noimport flag and that
#  paired sequences are in the panda directory, not the panda table in the database.
#  That means this script must be run with --load_seqs so the unique sequences are
#  loaded into the panda table.
#
#  TBD: make --load_seqs optional.  If not specified, it means the sequences to check
#  are in the panda table.  Write them in FASTA format in the working directory before
#  running cd-hit-dup.
#  ***************

import sqlite3
import argparse
import os
import os.path
import sys

from common import *
from FASTQ import *

###
# Run cd-hit-dup using options specified on the command line

input_file_pattern = 'merged.{}.fasta'
output_file_pattern = 'unique.{}.fasta'

def run_cd_hit_dup(sid, args):
    cmnd = 'cd-hit-dup'
    cmnd += ' -i ' + os.path.join(args.directory, input_file_pattern.format(sid))
    cmnd += ' -o ' + os.path.join(args.workspace, output_file_pattern.format(sid))
    print(cmnd)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
    db.commit()
    res = os.system(cmnd)
    
###
# Open the cluster file, get info about each cluster, return as a dictionary.
# The helper 'get_cluster' reads the lines in one cluster, returning the 
# line that follows a cluster and the lines in a cluster.  Save cluster info
# in the dictionary.

def get_cluster(line, file, info):
    cid = int(line.split()[1])        # first line is '>Cluster N'
    line = file.readline()
    seqinfo = line.split()[2]        # next line has defline of first seq in the cluster
    # defline = ':'.join(seqinfo[1:].split(':')[:-1])
    defline = FASTQ.parse_defline(seqinfo)
    count = 0
    while len(line) > 0 and line[0] != '>':
        count += 1
        line = file.readline()
    info[defline] = count
    return line

def parse_clusters(args, sid):
    clusters = { }
    file = open(os.path.join(args.workspace, output_file_pattern.format(sid) + '.clstr'))
    line = file.readline()
    while len(line) > 0:
        line = get_cluster(line, file, clusters)
    return clusters

###
# Load sequences from the fasta file produced by cd-hit-dup into the panda (paired
# sequence) table.

# Note: The panda table mirrors the uniq table in terms of which samples it holds.
# When this function is called the script is just about to load the unique sequences
# for a sample, so we need to delete any old sequences for that sample from the
# panda table.

insert_sequence = 'INSERT INTO panda (sample_id, defline, sequence) VALUES (?, ?, ?)'

def import_sequences(db, args, sid):
    db.execute('DELETE FROM panda WHERE sample_id = ?', (sid, ))
    file = open(os.path.join(args.workspace, output_file_pattern.format(sid)))
    defline = file.readline()
    while len(defline) > 0:
        seq = file.readline()
        defline = FASTQ.parse_defline(defline)
        db.execute(insert_sequence, (sid, defline, seq.strip()))
        defline = file.readline()

###
# Two things to do after each run of cd-hit-dup: if the panda table is empty we
# need to import the sequences, then import the clusters into the uniq table.

insert_cluster = 'INSERT INTO uniq (panda_id, sample_id, n) VALUES (?, ?, ?)'

def import_results(db, args, sid):
    if args.load_seqs:
        import_sequences(db, args, sid)
    defmap = defline_map(db, sid)
    clusters = parse_clusters(args, sid)
    for defline, size in clusters.items():
        db.execute(insert_cluster, (defmap[defline], sid, size))

###
# Top level function: run cd-hit-dup for each set of paired sequences, then save
# the id of the first sequence in the cluster and the number of sequences in the
# cluster.

def remove_duplicates(db, args):
    init_workspace(args)
    for row in sample_list(db, args):
        sid = row[0]
        run_cd_hit_dup(sid, args)
        import_results(db, args, sid)
            
###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    
    args = init_api(
        desc = "Run cd-hit-dup to find unique assembled sequences.",
        specs = [
            ('directory',    { 'metavar': 'dir', 'help' : 'name of directory containing assembled FASTQ files', 'default' : 'panda' } ),
            ('workspace',    { 'metavar': 'dir', 'help' : 'working directory', 'default' : 'uniq' } ),
            ('load_seqs',    { 'action': 'store_true', 'help' : 'load unique sequences into panda table'} ),
            ('sample',       { 'metavar': 'id', 'help' : 'process sequences from this sample only'} ),
        ]
    )
    
    if not args.load_seqs:          # see note at the front of this file
        print('This version requires --load_seqs (see documentation)')
        argparse.ArgumentParser.exit(1, 'Script aborted')
    
    db = sqlite3.connect(args.dbname)
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))

    try:
        uniq_spec = [('panda_id', 'foreign', 'panda'), ('sample_id', 'foreign', 'samples'), ('n', 'INTEGER')]
        init_table(db, 'uniq', 'uniq_id', uniq_spec, args.force, args.sample)
    except Exception as err:
        print('Error while initializing output table:', err)
        argparse.ArgumentParser.exit(1, 'Script aborted')

    remove_duplicates(db, args)
    record_metadata(db, 'end', '')
        
    db.commit()
