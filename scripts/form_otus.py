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

# Define filenames (may be referenced by more than one function)

input_file = 'seeds.fasta'
cluster_file = 'clusters.fasta'

###
# Print the sequences that will be clustered.  If a prior stage did a map to reference
# sequences the results will be in a table named 'hits' and we need to merge that table
# with the unique sequences.  Otherwise just print the unique sequences in order of 
# decreasing frequency
    
find_ref_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "hits"'

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
    record_metadata(db, 'query', sql)
    return db.execute(sql)
    
def print_unique_sequences(db, args):
    res = fetch_unique_sequences(db, args)
    ff = open(os.path.join(args.workspace, input_file), 'w')
    for pid, n, defline, sequence in res:
        print_one_seq(ff, n, defline, sequence)
    ff.close()

select_hits = 'SELECT panda_id, identity, match_id, match_chars FROM hits'

def merge_ref_with_unique(db, args):
    record_metadata(db, 'query', select_hits)
    
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
    record_metadata(db, 'exec', cmnd, commit=True)
    res = os.system(cmnd)
    
###
# Populate the tables by importing the FASTA file produced by uclust.  

insert_cluster = 'INSERT INTO clusters (cluster_id, name, sequence) VALUES (?,?,?)'
insert_member = 'INSERT INTO members (cluster_id, name, diffsqm, dqt, dqm) VALUES (?,?,?,?,?)'

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
    
if __name__ == "__main__":
    
    args = init_api(
        desc = "Run uclust to create de novo OTUs.",
        specs = [
            ('workspace',    { 'metavar': 'dir', 'help' : 'working directory', 'default' : 'clusters' } ),
            ('singletons',   { 'action': 'store_true', 'help' : 'include singletons (default: disregard singletons)' } ),
        ]
    )
        
    db = sqlite3.connect(args.dbname)
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))

    try:
        cluster_spec = [('name', 'TEXT'), ('sequence', 'TEXT')]
        init_table(db, 'clusters', 'cluster_id', cluster_spec, args.force)
        member_spec = [('cluster_id', 'foreign', 'clusters'), ('name', 'TEXT'), ('diffsqm', 'INTEGER'), ('dqt', 'REAL'), ('dqm', 'REAL')]
        init_table(db, 'members', 'member_id', member_spec, args.force)
    except Exception as err:
        print('Error while initializing output tables:', err)
        argparse.ArgumentParser.exit(1, 'Script aborted')

    form_otus(db, args)
    record_metadata(db, 'end', '')

    db.commit()
