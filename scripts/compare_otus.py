#! /usr/bin/env python3

# Compare clusters created by two different 16S rRNA pipelines

# John Conery
# University of Oregon
# 2014-12-08

# Usage:
#
#    compare_otus.py [--info] DB1 DB2
#

import argparse
import sqlite3
import os.path
import sys
import re

##
# Execute a query that returns a single value, e.g. a table size, or None if the result
# set is empty

def fetch_value(db, query):
    res = db.execute(query).fetchall()
    if len(res) > 0:
        return res[0][0]
    else:
        return None
        
##
# Return the number of items in the clusters tables

def cluster_table_sizes(db, x):
    
    def table_size(t):
        sql = 'SELECT count(*) FROM {dbname}.{tbl}'.format(dbname=x, tbl=t)
        # print(sql)
        return fetch_value(db, sql)
        
    return list(map(table_size,['clusters', 'members']))
    
##
# Print output in tabular form, computing column sizes automatically

def col_size(res):
    maxlen = lambda i: max(map(len,map(lambda tup: str(tup[i]), res)))
    widths = map(maxlen,range(len(res[0])))
    linespec = "  ".join(map(lambda n: "{:%d}" % n, widths))
    return linespec
    
def print_header(names, linespec):
    seps = map(lambda x: '-'*len(x), names)
    print(linespec.format(*names))
    print(linespec.format(*seps))
    
## 
# Prep:  make tables showing clusters in both DBs and in one DB but not the other

a_and_b = 'CREATE TABLE AandB AS SELECT A.clusters.cluster_id as aid, A.clusters.name AS aname, B.clusters.cluster_id as bid, B.clusters.name AS bname FROM A.clusters JOIN B.clusters USING (sequence)'
a_not_b = 'CREATE TABLE AnotB AS SELECT A.clusters.cluster_id as aid, A.clusters.name AS aname FROM A.clusters LEFT JOIN B.clusters USING (sequence) WHERE b.clusters.name IS NULL'
b_not_a = 'CREATE TABLE BnotA AS SELECT B.clusters.cluster_id as bid, B.clusters.name AS bname FROM B.clusters LEFT JOIN A.clusters USING (sequence) WHERE a.clusters.name IS NULL'

taxa = 'CREATE TABLE {d}taxa AS SELECT cluster_id, genus.name AS genus FROM {d}.clusters LEFT JOIN {d}.taxonomy ON (cluster_id = otu_id) JOIN {d}.genus USING (genus_id)'

def make_temp_tables(db):
    db.execute(a_and_b)
    db.execute(a_not_b)
    db.execute(b_not_a)
    
    db.execute(taxa.format(d='A'))
    db.execute(taxa.format(d='B'))
    
##
# Print summary info about clusters in each database

def cluster_info(db, args):
    linespec = '{:15s} {:10d} {:10d}'
    print('                      clusters    members')
    print('                      --------    -------')
    print(linespec.format(args.db1, *cluster_table_sizes(db, 'A')))
    print(linespec.format(args.db2, *cluster_table_sizes(db, 'B')))
    print()
    
    print('Clusters in both {} and {}: '.format(args.db1, args.db2), fetch_value(db, 'SELECT count(*) FROM AandB'))
    print('Clusters in {} but not {}:  '.format(args.db1, args.db2), fetch_value(db, 'SELECT count(*) FROM AnotB'))
    print('Clusters in {} but not {}:  '.format(args.db2, args.db1), fetch_value(db, 'SELECT count(*) FROM BnotA'))

##
# Print more detailed info about each cluster

def cluster_details(db):
    for aid, bid in db.execute('SELECT aid, bid FROM AandB'):
        print(aid, 'vs.', bid, file=sys.stderr)

        # sql = 'SELECT sequence from AandB join A.members on (aid = A.members.cluster_id) join A.panda on (defline = name) where aid = {}'.format(aid)
        # res = db.execute(sql)
        # seqs_in_a = set(map(lambda x: x[0], res))
        #
        # sql = 'SELECT sequence from AandB join B.members on (bid = B.members.cluster_id) join B.panda on (defline = name) where bid = {}'.format(bid)
        # res = db.execute(sql)
        # seqs_in_b = set(map(lambda x: x[0], res))
        
        seqs_in_a = sequence_set(db, 'A', aid)
        seqs_in_b = sequence_set(db, 'B', bid)
        
        ataxon = fetch_value(db, 'SELECT genus FROM Ataxa WHERE cluster_id = {}'.format(aid))
        btaxon = fetch_value(db, 'SELECT genus FROM Btaxa WHERE cluster_id = {}'.format(bid))
        
        # print('\t'.join(map(str,[aid, bid, ataxon, btaxon, len(seqs_in_a & seqs_in_b), len(seqs_in_a - seqs_in_b), len(seqs_in_b - seqs_in_a)])))
        print_line(aid, bid, ataxon, btaxon, seqs_in_a, seqs_in_b)

    for aid, aname in db.execute('SELECT aid, aname FROM AnotB'):
        seqs_in_a = sequence_set(db, 'A', aid)
        ataxon = fetch_value(db, 'SELECT genus FROM Ataxa WHERE cluster_id = {}'.format(aid))
        print_line(aid, None, ataxon, None, seqs_in_a, set())
        
    for bid, bname in db.execute('SELECT bid, bname FROM BnotA'):
        seqs_in_b = sequence_set(db, 'B', bid)
        btaxon = fetch_value(db, 'SELECT genus FROM Btaxa WHERE cluster_id = {}'.format(bid))
        print_line(None, bid, None, btaxon, set(), seqs_in_b)

def sequence_set(db, x, row):
    sql = 'SELECT sequence from {d}.members join {d}.panda on (defline = name) where cluster_id = {r}'.format(d=x, r=row)
    res = db.execute(sql)
    return set(map(lambda x: x[0], res))

def print_line(aid, bid, atax, btax, aset, bset):
    print('\t'.join(map(str,[aid, bid, atax, btax, len(aset & bset), len(aset - bset), len(bset - aset)])))
    

###
# Parse command line arguments...

def init_api():
    parser = argparse.ArgumentParser(
        description="...",
        epilog="..."
    )
    parser.add_argument('db1', help='baseline database name')
    parser.add_argument('db2', help='comparison database name')
    parser.add_argument('-i', '--info', help='print cluster sizes and exit', action='store_true')
    return parser.parse_args()

###
# Top level....
    
if __name__ == "__main__":
    args = init_api()
    db = sqlite3.connect(':memory:')
    db.execute('attach database "{}" as A'.format(args.db1))
    db.execute('attach database "{}" as B'.format(args.db2))
    make_temp_tables(db)
    if args.info:
        cluster_info(db, args)
    else:
        cluster_details(db)
