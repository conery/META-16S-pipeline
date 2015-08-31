#! /usr/bin/env python3

# Run PANDAseq to find overlapping parts of the  two ends of a paired-end read.  PANDAseq
# also does spacer and primer trimming and quality filtering.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-29


import sqlite3
import argparse
import os
import os.path
import sys

from common import *
from FASTQ import *

# File names used in calls to pandaseq

log_file_pattern = 'log.{}.txt'
merge_file_pattern = 'merged.{}.fasta'

###
# Fetch the primer sequences from the database

def fetch_primers(db):
    """
    Read the primer sequences from the database
    """
    primers = { }
    for n, seq in db.execute("SELECT pair_end, sequence FROM primers"):
        primers[int(n)] = seq
    return primers

###
# Run PANDAseq using options specified on the command line
# TBD: write to /dev/null?  -N??

def run_panda(args, sid, fn1, fn2, primers):
    cmnd = 'pandaseq'
    if args.algorithm:
        cmnd += ' -A ' + args.algorithm
    cmnd += ' -N'                           # throw out seqs with N's
    cmnd += ' -f ' + os.path.join(args.directory, fn1)
    cmnd += ' -r ' + os.path.join(args.directory, fn2)
    if len(primers) > 0:
        cmnd += ' -p ' + primers[1]
        cmnd += ' -q ' + primers[2]
    cmnd += ' -w ' + os.path.join(args.workspace, merge_file_pattern.format(sid))
    cmnd += ' -g ' + os.path.join(args.workspace, log_file_pattern.format(sid))
    if args.minlength is not None:
        cmnd += ' -l ' + str(args.minlength)
    if args.minoverlap is not None:
        cmnd += ' -o ' + str(args.minoverlap)
    print(cmnd)
    if not args.norun:
        record_metadata(db, 'exec', cmnd, commit=True)
        res = os.system(cmnd)

###
# Import the assembled sequences

insert_sequence = 'INSERT INTO panda (sample_id, defline, sequence) VALUES (?, ?, ?)'

def import_results(db, args, sid):
    file = open(os.path.join(args.workspace, merge_file_pattern.format(sid)))
    defline = file.readline()
    count = 0
    limit = int(args.limit) if args.limit else None
    while len(defline) > 0:
        defline = FASTQ.parse_defline(defline)
        sequence = file.readline()
        db.execute(insert_sequence, (sid, defline, sequence.strip()))
        count += 1
        if args.limit and count >= limit:
            break
        defline = file.readline()

###
# Top level function: initialize the workspace directory, get sample parameters from
# the database, run PANDAseq for all specified samples

def assemble_pairs(db, args):
    init_workspace(args)
    primers = fetch_primers(db)
    for sid, sname, fn1, fn2 in sample_list(db, args):
        run_panda(args, sid, fn1, fn2, primers)
        if not args.noimport and not args.norun:
            import_results(db, args, sid)

###
# Check the combination of command line options to make sure they're sensible

def validate_options(args):
    # if we're not loading data the limit option is superfluous
    if args.noimport and args.limit:
        print('Warning: options ignored: with --noimport the following options are ignored: --limit')

###
# Parse the command line arguments, call the top level function...

if __name__ == "__main__":
    
    args = init_api(
        desc = "Run PANDAseq to filter low quality sequences, and find overlapping ends of pairs of reads. The default algorithm is pandaseq; alternatives are simple_bayesian, ea_util, flash, pear, rdp_mle, stitch, and uparse.",
        with_limits = True,
        specs = [
            ('directory',    { 'metavar': 'dir', 'help' : 'name of directory containing FASTQ files' } ),
            ('workspace',    { 'metavar': 'dir', 'help' : 'working directory', 'default' : 'panda' } ),
            ('algorithm',    { 'metavar': 'name', 'help' : 'algorithm for computing overlaps'} ),
            ('minlength',    { 'metavar': 'N', 'help' : 'minimum assembled sequence length', 'type' : int} ),
            ('minoverlap',   { 'metavar': 'N', 'help' : 'minimum overlap length', 'type' : int} ),
            ('norun',        { 'action': 'store_true', 'help' : "print shell commands but don't execute them"} ),
        ]
    )
    
    validate_options(args)
    db = sqlite3.connect(args.dbname)
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))
    
    try:
        panda_spec = [('sample_id', 'foreign', 'samples'), ('defline', 'TEXT'),  ('sequence', 'TEXT')]
        init_table(db, 'panda', 'panda_id', panda_spec, args.force, args.sample)
    except Exception as err:
        print('Error while initializing output table:', err)
        argparse.ArgumentParser.exit(1, 'Script aborted')
    
    assemble_pairs(db, args)
    record_metadata(db, 'end', '')
    
    db.commit()
