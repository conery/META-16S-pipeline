#! /usr/bin/env python3

# Import reads into a SQLite database for a PIP/NGS analysis pipeline.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-15

# Read a set of fastq files, extract the index from the defline, pack
# the sequence and quality chars into a "BLOB"

import sqlite3
import argparse
import os.path
import sys

from common import *
from FASTQ import *

###
# Read the sequences from a pair of files (R1 and R2) for a specified sample.  Create
# a FASTQ object for each read, compress it, and insert the compressed form into the 
# reads table.

insert_blob = 'INSERT INTO reads (sample_id, read) VALUES (?, ?)'
insert_blob_and_def = 'INSERT INTO reads (sample_id, read, defline) VALUES (?, ?, ?)'
    
def load_sequences(db, fn1, fn2, sid, args):
    "Load files fn1 and fn2 into the reads table"
    
    if not (os.path.exists(fn1) and os.path.exists(fn2)):
        return
    
    file1 = FASTQReader(fn1)
    file2 = FASTQReader(fn2)
    count = 0
    
    limit = int(args.limit) if args.limit else None

    seq1 = file1.readseq()
    seq2 = file2.readseq()
    
    while seq1 is not None and seq2 is not None:
        for seq in [seq1, seq2]:
            if args.quality and seq.filtered() == 'Y':
                blob = None
            else:
                seq.pack()
                blob = seq.blob()
            if args.deflines:
                # db.execute(insert_blob_and_def, (sid, blob, seq.unique_defline()))
                db.execute(insert_blob_and_def, (sid, blob, seq.defline()))
            else:
                db.execute(insert_blob, (sid, blob))
        count += 1
        if limit and count >= limit:
            break
        seq1 = file1.readseq()
        seq2 = file2.readseq()
        
    file1.close()
    file2.close()
    
###
# The main function -- iterate over the names of the samples, call
# load_sequences to import the sequences from the FASTQ files.
    
def import_files(db, samples, args):
    for sid, sname, fastq1, fastq2 in samples:
        fn1 = os.path.join(args.directory, fastq1)
        fn2 = os.path.join(args.directory, fastq2)
        if not args.noimport:
            load_sequences(db, fn1, fn2, sid, args)
            record_metadata(db, 'import', '{}, {}'.format(fn1,fn2))

    # TBD: consider making the index a command line option
    if not args.noimport:
        db.execute('CREATE INDEX IF NOT EXISTS defx ON reads (defline)')

###
# Check the combination of command line options to make sure they're sensible

def validate_options(args):
    # if we're not loading data the limit, defline, and quality options are superfluous
    if args.noimport and (args.limit or args.deflines or args.quality):
        print('Warning: options ignored: with --noimport the following options are ignored: --limit, --defline, --quality')

###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    
    # Make the set of arguments common to all apps (including --limit, --sample, and --noimport), then add
    # args specific to this script
    args = init_api(
        desc = "Import FASTQ files into a SQLite3 database for a PIP/NGS analysis pipeline.",
        epi = "If --noimport is specified sequence descriptions are created but sequences are not loaded into the database.",
        with_limits = True,
        specs = [
            ('directory',    { 'required' : True, 'metavar': 'dir', 'help' : '(required) name of directory containing FASTQ files' } ),
            ('quality', { 'action': 'store_true', 'help' : "don't import sequences flagged as low quality"} ),
            ('deflines', { 'action': 'store_true', 'help' : 'include deflines with sequences'} ),
        ]
    )
    
    validate_options(args)
    db = sqlite3.connect(args.dbname)
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))
    
    try:
        read_spec = [('sample_id', 'foreign', 'samples'), ('read', 'BLOB'),  ('defline', 'TEXT DEFAULT NULL')]
        init_table(db, 'reads', 'read_id', read_spec, args.force, args.sample)
    except Exception as err:
        print('Error while initializing output table:', err)
        argparse.ArgumentParser.exit(1, 'Script aborted')
    
    samples = sample_list(db, args)                
    import_files(db, samples, args)
        
    record_metadata(db, 'end', '')
    db.commit()
