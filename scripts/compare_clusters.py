#! /usr/bin/env python3

# Compare clusters created by two different 16S rRNA pipelines.  Uses vsearch to align
# the central sequences of clusters in one database with those in a second database.

# John Conery
# University of Oregon
# 2015-08-25

# Usage:
#
#    compare_clusters.py DB1 DB2
#

import argparse
import sqlite3
import os.path
import sys
import re
import tempfile
import subprocess

##
# Execute a query that returns a single value, e.g. a table size, or None if the result
# set is empty

def fetch_value(db, query):
    res = db.execute(query).fetchall()
    if len(res) > 0:
        return res[0][0]
    else:
        return None

###
# Get the number of sequences in each cluster table

def get_counts(db):
    res = []
    for tbl in ['A.clusters', 'B.clusters']:
        sql = 'SELECT count(*) FROM {}'.format(tbl)
        if args.minlength:
            sql += ' WHERE length(sequence) > {}'.format(args.minlength)
        res.append(fetch_value(db,sql))
    return res

###
# Write sequences to a FASTA file, return the name of the file.

def write_fasta(tbl, args):
    fn = tbl + '.fa'
    with open(fn, 'w') as f:
        sql = 'SELECT name, sequence FROM {}'.format(tbl)
        if args.minlength is not None:
            sql += ' WHERE length(sequence) > {}'.format(args.minlength)
        # print(sql)
        for res in db.execute(sql).fetchall():
            print('>', res[0], file=f)
            print(res[1], file=f)
    return fn
 
###
# Run vsearch

search_cmnd = 'vsearch --usearch_global {} -db {} -id 0.97 -strand plus -userout {} -userfields query+target+id+ql+pairs'

def run_vsearch(f1, f2, out):
    cmnd = search_cmnd.format(f1, f2, out)
    # print(cmnd)
    subprocess.call(cmnd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

###
# The number of lines in an output file is the number of sequences found in the other database.

def number_of_hits(fn):
    res = subprocess.check_output('wc {}'.format(fn), shell=True)
    return int(res.strip().split()[0])

###
# Compare each file against the other.  

def compare_with_vsearch(args, na, nb):
    fa1 = write_fasta('A.clusters', args)
    fa2 = write_fasta('B.clusters', args)
    try:
        run_vsearch(fa1, fa2, 'ab.tsv')
        ha = number_of_hits('ab.tsv')
        run_vsearch(fa2, fa1, 'ba.tsv')
        hb = number_of_hits('ba.tsv')
    except Exception as err:
        print(err)
    print(args.db1, '=>', args.db2+':', ha, '/', na)
    print(args.db2, '=>', args.db1+':', hb, '/', nb)
    print('similarity: %5.2f' % ((ha + hb)/(na + nb)))
    

###
# Parse command line arguments...

def init_api():
    parser = argparse.ArgumentParser(
        description="Compare 'clusters' tables in two databases.",
    )
    parser.add_argument('db1', help='baseline database name')
    parser.add_argument('db2', help='comparison database name')
    parser.add_argument('-m', '--minlength', help='minimum sequence length', type=int)
    parser.add_argument('-i', '--info', help='print table sizes and exit', action='store_true')
    return parser.parse_args()

###
# Top level....
    
if __name__ == "__main__":
    args = init_api()
    db = sqlite3.connect(':memory:')
    db.execute('attach database "{}" as A'.format(args.db1))
    db.execute('attach database "{}" as B'.format(args.db2))
    na, nb = get_counts(db)
    if args.info:
        print('clusters in {}: '.format(args.db1), na)
        print('clusters in {}: '.format(args.db2), nb)
    else:
        compare_with_vsearch(args, na, nb)
