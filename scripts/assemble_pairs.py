#! /usr/bin/env python3

# Run PANDAseq to find overlapping parts of the  two ends of a paired-end read.  PANDAseq
# also does spacer and primer trimming and quality filtering.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-29

from FASTQ import *

import sqlite3
import argparse
import os
import os.path
import sys

# These global vars define filenames and other items shared by two or more steps

table_name = 'panda'
log_file_pattern = 'log.{}.txt'
merge_file_pattern = 'merged.{}.fasta'

###
# Initialize (or re-initialize) the table that will hold the output from PANDAseq

find_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{}"'.format(table_name)
drop_table = 'DROP TABLE {}'.format(table_name)
create_table = 'CREATE TABLE {tbl} ( {tbl}_id INTEGER PRIMARY KEY AUTOINCREMENT, sample_id INTEGER, defline TEXT, sequence TEXT )'.format(tbl = table_name)

def prepare_table(db, args):
    """
    Return True if the 'panda' table is initialized and ready to accept values.  If the
    table exists already don't overwrite it unless --force was specified on the command line.
    """
    if db.execute(find_table).fetchall():
        if args.force:
            db.execute(drop_table)
        else:
            return False
    db.execute(create_table)
    return True

###
# Initialize the directory where intermediate work products will be stored.

def init_workspace(args):
    dirname = args.workspace
    
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    else:
        for root, dirs, files in os.walk(dirname, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))

###
# Fetch sample names and the corresponding file names

# fetch_barcodes = 'SELECT barcode_id, experiment, sequence FROM barcodes'
# # file_pattern = 'L1_{exp}_{code1}-{code2}_L001_R{pair}_001.fastq'
#
# def load_barcodes(db, args):
#     "Generate the names of the FASTQ files, associate them with barcode IDs."
#     files = { }
#     sql = fetch_barcodes
#     if isinstance(args.experiment, str):
#         sql += ' WHERE experiment = "{}"'.format(args.experiment)
#     for bcid, expt, code in db.execute(sql):
#         fn1 = os.path.join(args.directory, file_pattern.format(exp = expt, code1 = code[:7], code2 = code[7:], pair = 1))
#         fn2 = os.path.join(args.directory, file_pattern.format(exp = expt, code1 = code[:7], code2 = code[7:], pair = 2))
#         files[bcid] = {1: fn1, 2: fn2}
#     return files
    
fetch_samples = 'SELECT sample_id, name, r1_file, r2_file FROM samples'

def fetch_sample_list(db, args):
    "Return a list of sample names and the corresponding FASTQ file names."
    sql = fetch_samples
    if isinstance(args.sample, str):
        sql += ' WHERE name = "{}"'.format(args.sample)
    return db.execute(sql).fetchall()

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
    # cmnd += ' -d B'                       # debug/trace info to keep
    cmnd += ' -N'                           # throw out seqs with N's
    # cmnd += ' -t 0.7'                       # quality threshold (default = 0.6)
    cmnd += ' -f ' + os.path.join(args.directory, fn1)
    cmnd += ' -r ' + os.path.join(args.directory, fn2)
    if len(primers) > 0:
        cmnd += ' -p ' + primers[1]
        cmnd += ' -q ' + primers[2]
    # else:
    #     cmnd += ' -A flash'
    cmnd += ' -w ' + os.path.join(args.workspace, merge_file_pattern.format(sid))
    cmnd += ' -g ' + os.path.join(args.workspace, log_file_pattern.format(sid))
    # cmnd += ' -l ' + str(args.minlength)
    cmnd += ' -o ' + str(args.minoverlap)
    print(cmnd)
    if not args.norun:
        db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
        db.commit()
        res = os.system(cmnd)

###
# Import the assembled sequences

insert_sequence = 'INSERT INTO {} (sample_id, defline, sequence) VALUES (?, ?, ?)'.format(table_name)

def import_results(db, args, sid):
    file = open(os.path.join(args.workspace, merge_file_pattern.format(sid)))
    defline = file.readline()
    count = 0
    limit = int(args.limit) if args.limit else None
    while len(defline) > 0:
        # defline = ':'.join(defline[1:].split(':')[:-1])
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
    # files = fetch_fastq_filenames(db, args)
    for sid, sname, fn1, fn2 in fetch_sample_list(db, args):
        run_panda(args, sid, fn1, fn2, primers)
        if not args.noimport:
            import_results(db, args, sid)

###
# Set up command line arguments

def init_api():
    parser = argparse.ArgumentParser(
        description="""Run PANDAseq to remove spacers and primers, filter low quality sequences, and find
        the overlapping ends of pairs of reads.
        """
    )
    parser.add_argument('dbname', help='the name of the SQLite database file')
    parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
    parser.add_argument('-w', '--workspace', help='working directory', default='panda')
    parser.add_argument('-d', '--directory', help='directory with input FASTQ files', default='data')
    parser.add_argument('-s', '--sample', metavar='id', required=False, help='sample to assemble')
    parser.add_argument('-a', '--all', action='store_true', help='assemble all samples')
    # parser.add_argument('-m', '--minlength', help='minimum assembled sequence length', type=int, default=150)
    parser.add_argument('-m', '--minoverlap', help='minimum overlap length', type=int, default=10)
    parser.add_argument('-L', '--limit', help='number of sequences to import')
    parser.add_argument('--noimport', action='store_true', help="assemble pairs, but leave them in folder")
    parser.add_argument('--norun', action='store_true', help="print shell commands but don't execute them")
    return parser.parse_args()

###
# Parse the command line arguments, call the top level function...

if __name__ == "__main__":
    args = init_api()
    
    if not (args.sample or args.all):
        argparse.ArgumentParser.exit(1, 'Specify a sample name or --all for all samples')
        
    if args.limit and args.noimport:
        argparse.ArgumentParser.exit(1, 'Incompatible options:  --limit and --noimport')
    
    db = sqlite3.connect(args.dbname)
    
    if not prepare_table(db, args):
        argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')
    
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
    assemble_pairs(db, args)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )
    
    db.commit()
