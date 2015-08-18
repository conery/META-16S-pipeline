#! /usr/bin/env python3

# Run usearch to create de novo OTUs.  Each OTU will be represented by the 
# "seed" sequence used by uclust to as the center of each OTU.

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-04-15

from common import *

import sqlite3
import argparse
import os
import os.path
import re
import sys

# These global vars define filenames and other items shared by two or more steps

ref_table_name = 'hits'
cluster_table_name = 'clusters'
member_table_name = 'members'
input_file = 'seeds.fasta'
cluster_file = 'clusters.fasta'
# otu_file = 'otus.fasta'

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
# Initialize (or re-initialize) the tables that will hold the output from uclust

find_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{tbl}"'
drop_table = 'DROP TABLE {tbl}'
create_table = {
    cluster_table_name : 'CREATE TABLE {tbl} ( cluster_id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, sequence TEXT )'.format(tbl = cluster_table_name),
    member_table_name : 'CREATE TABLE {tbl} ( member_id INTEGER PRIMARY KEY AUTOINCREMENT, cluster_id INTEGER, name TEXT, diffsqm INTEGER, dqt REAL, dqm REAL )'.format(tbl=member_table_name),
}

def prepare_tables(db, args):
    """
    Return True if the results tables are initialized and ready to accept values.  If a
    table exists already don't overwrite it unless --force was specified on the command line.
    """
    for tbl in [cluster_table_name, member_table_name]:
        if db.execute(find_table.format(tbl=tbl)).fetchall():
            if args.force:
                db.execute(drop_table.format(tbl=tbl))
            else:
                return False
        db.execute(create_table[tbl])
    return True

###
# Print sequences in the uniq table in fasta format

# Update:  first see if there is a 'hits' table with results of the map-to-reference step.
# If so, select best hits, use them to create the first part of the seeds file; if not,
# just initialize an empty seeds file. The code that writes the unique sequecnes will 
# append them to the seeds file.


# def init_seed_file(db, args):
#         db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'query', make_best_hits))
#         db.execute(make_best_hits)
#         db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'query', fetch_hits))
#         ff = open(os.path.join(args.workspace, input_file), 'w')
#         for query, target, n, pct, seq in db.execute(fetch_hits).fetchall():
#             print('>{};size={}'.format(target,n), file=ff)
#             print(re.sub('-','',seq), file=ff)
#         ff.close()
#     else:
#         ff = open(os.path.join(args.workspace, input_file), 'w')
#         ff.close()
        
###
# Print the sequences that will be clustered.  If a prior stage did a map to reference
# sequences the results will be in a table named 'hits' and we need to merge that table
# with the unique sequences.  Otherwise just print the unique sequences in order of 
# decreasing frequency
    
find_ref_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{}"'.format(ref_table_name)

def print_sequences(db, args):
    if db.execute(find_ref_table).fetchall():
        merge_ref_with_unique(db, args)
    else:
        print_unique_sequences(db, args)
        
# select_unique_sequences = 'SELECT panda_id, n, defline, sequence FROM uniq JOIN panda USING (panda_id)'

def fetch_unique_sequences(db, args):
    sql = 'SELECT panda_id, sum(n) as count, defline, sequence FROM uniq JOIN panda USING (panda_id)'
    if not args.singletons:
        sql += ' WHERE n > 1'
    sql += ' GROUP by sequence'
    sql += ' ORDER by count DESC'
    # print(sql)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'query', sql))
    return db.execute(sql)
    
def print_unique_sequences(db, args):
    res = fetch_unique_sequences(db, args)
    ff = open(os.path.join(args.workspace, input_file), 'w')
    for pid, n, defline, sequence in res:
        print_one_seq(ff, n, defline, sequence)
        # print('>{};size={}'.format(defline,n), file=ff)
        # print(sequence, file=ff)
    ff.close()

# Deprecated -- hits table from map_reference.py now has only one best hit
# make_index = 'CREATE INDEX qx ON hits (query)'
# make_best_hits = 'CREATE TEMPORARY TABLE best_hits AS SELECT query, target, max(identity) AS identity, CAST(align_length as real) / query_length AS align_pct, target_chars FROM hits GROUP BY query'
# select_hits = 'SELECT n, query, target_chars FROM panda JOIN uniq USING (panda_id) JOIN best_hits ON (query = defline) WHERE align_pct > 0.95 ORDER BY n DESC'

select_hits = 'SELECT panda_id, identity, match_id, match_chars FROM {tbl}'.format(tbl=ref_table_name)

def merge_ref_with_unique(db, args):
    # print(select_hits)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'query', select_hits))
    
    hits = { }
    for pid, pct, match_id, match_chars in db.execute(select_hits):
        hits[pid] = { 'ident' : pct, 'id' : match_id, 'chars' : match_chars}
    
    ff = open(os.path.join(args.workspace, input_file), 'w')
    
    for pid, n, defline, sequence in fetch_unique_sequences(db, args):
        if pid in hits:
            target = hits[pid]
            if target['ident'] < 100.0:
                print_one_seq(ff, n, target['id'], re.sub('-','',target['chars']), ref=True)
        print_one_seq(ff, n, defline, sequence)
        
    ff.close()

def print_one_seq(f, n, defline, seq, ref=False):
    if ref:
        defline += ':ref'
    print('>{};size={}'.format(defline,n), file=f)
    print(seq, file=f)    
    
###
# Run the app.  

def run_cluster_otus(args):
    cmnd = 'usearch -cluster_otus '
    cmnd += os.path.join(args.workspace, input_file)
    # cmnd += ' -otus ' + os.path.join(args.workspace, otu_file)
    cmnd += ' -fastaout ' + os.path.join(args.workspace, cluster_file)
    print(cmnd)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
    db.commit()
    res = os.system(cmnd)
    
###
# Populate the tables by importing the FASTA file produced by uclust.  

insert_cluster = 'INSERT INTO {} (cluster_id, name, sequence) VALUES (?,?,?)'.format(cluster_table_name)
insert_member = 'INSERT INTO {} (cluster_id, name, diffsqm, dqt, dqm) VALUES (?,?,?,?,?)'.format(member_table_name)

def import_results(db, args):
    f = open(os.path.join(args.workspace, cluster_file))
    defline = f.readline()
    cid = 0
    cmap = { }

    while defline:
        sequence = ''
        seqline = f.readline()
        while seqline:
            if seqline.startswith('>'):  break
            sequence += seqline.strip()
            seqline = f.readline()
        # m = re.search(r'>(.*);size=(\d+)', defline)
        # db.execute(insert_otu, (dm[m.group(1)], sequence))
        d = parse_defline(defline.strip('>;\n'))
        if d['up'] == 'otu':
            db.execute(insert_cluster, (cid, d['name'], sequence))
            cmap[d['name']] = cid
            cid += 1
        elif d['up'] == 'member' and ':ref' not in d['name']:
            db.execute(insert_member, (cmap[d['top']], d['name'], d['diffsqm'], d['dqt'], d['dqm']))
        defline = seqline

def parse_defline(s):
    parts = s.split(';')
    d = {'name' : parts[0]}                # first item is sequence name, remaining are k=v pairs
    d.update(dict(map(lambda x: x.split('='), parts[1:])))
    return d
    
###
# Top level function: initialize the workspace directory, run the app

def form_otus(db, args):
    init_workspace(args)
    print_sequences(db, args)
    run_cluster_otus(args)
    import_results(db, args)

###
# Parse the command line arguments, call the top level function...

def init_api():
    parser = argparse.ArgumentParser(
        description="Run uclust to create de novo OTUs.",
    )
    parser.add_argument('dbname', help='the name of the SQLite database file')
    parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
    parser.add_argument('--singletons', action='store_true', help = 'include singletons (default: disregard singletons)')
    parser.add_argument('--workspace', help='working directory', default='clusters')
    return parser.parse_args()
    
###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    args = init_api()

    db = sqlite3.connect(args.dbname)
        
    if not prepare_tables(db, args):
        argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')

    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
    form_otus(db, args)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )

    db.commit()
