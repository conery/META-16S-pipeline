#! /usr/bin/env python3

# Determine the taxonomic classification of a set of sequences

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-06-05

#  **** NOTE ****
#  This version assumes the OTUs produced by the map_otus script are still in a
#  FASTA files in the working directory of that script.  The default directory 
#  name is 'otus' but an alternative can be specified with --directory.
#
#  TBD: don't have a default directory; if no directory is specified write the 
#  OTU sequences in the otus table to the working directory for this script.
#  ***************

import sqlite3
import argparse
import os
import os.path
from string import punctuation
import sys

from common import *

# These strings are the names of the columns in the taxononmy table.  Note that
# 'orderx' has an 'x' at the end -- 'order' is a keyword in SQL so we can't use
# it for a column name.

levels = ['domain', 'phylum', 'class', 'orderx', 'family', 'genus']

# These global vars define filenames and other items shared by two or more steps

input_file = 'otus.fasta'
output_file = 'classifications.txt'

###
# Run the app

def run_classifier(args):
    cmnd = 'java -jar ' + classifier_path
    cmnd += ' classify ' + os.path.join(args.directory, input_file)
    cmnd += ' -f fixrank'
    cmnd += ' -o ' + os.path.join(args.workspace, output_file)
    print(cmnd)
    record_metadata(db, 'exec', cmnd, commit=True)
    res = os.system(cmnd)

###
# Populate the table.  The function uses separate maps for each taxonomic level.  As new
# names are discovered they are added to the map with a unique integer key.  After
# reading all the assignments the maps are written as new tables in the database.

def import_results(db, args):
    
    names = dict(zip(levels, [dict() for i in range(len(levels))]))

    def sanitize(recs):
        'change "order" to "orderx", strip quotes'
        for i in range(len(recs)):
            if recs[i] == 'order':
                recs[i] = 'orderx'
            recs[i] = recs[i].strip(punctuation)

    def update_names(parts):
        'iterate over the parts of a taxonomy string to update name tables'
        for p, x in parts.items():          # p is a level (domain, phylum, etc), x is the name at that level
            if x is None or len(x) == 0: 
                continue
            d = names[p]
            if x not in d:
                d[x] = len(d)+1

    for line in open(os.path.join(args.workspace, output_file)):
        parts = dict(zip(levels, [None]*len(levels)))
        scores = dict(zip(levels, [0]*len(levels)))
        recs = line.split('\t')
        sanitize(recs)
        for i in range(2, len(recs), 3):
            parts[recs[i+1]] = recs[i]
            scores[recs[i+1]] = recs[i+2]
        update_names(parts)
        save_record(recs[0], names, parts, scores)

    for x in levels:
        db.execute('DROP TABLE IF EXISTS {}'.format(x))
        db.execute('CREATE TABLE {} ( {}_id INTEGER PRIMARY KEY, name TEXT)'.format(x,x))
    for x, values in names.items():
        for name, item_id in values.items():
            q = 'INSERT INTO {} ({}_id, name) VALUES (?,?)'.format(x,x)
            db.execute(q, (item_id, name))
    

###
# Helper function for import_results -- write a record to the taxonomy table using the
# values in the names dictionry as foreign keys

cols = ', '.join(map(lambda x: x + '_id', levels))
cols += ', ' + ', '.join(map(lambda x: 'p_' + x, levels))
fmts = ','.join(map(lambda x: '{' + x +'}', levels))
pfmts = ','.join(map(lambda x: '{p' + x +'}', levels))
insert_record = 'INSERT INTO taxonomy (otu_id, {}) VALUES (?,{},{})'.format(cols, fmts, pfmts)

def save_record(otu, counts, parts, scores):
    d = {}
    c = {}
    for x in parts:
        if parts[x] is None or parts[x] == '':
            d[x] = 'NULL'
            c['p'+x] = 'NULL'
        else:
            d[x] = counts[x][parts[x]]
            c['p'+x] = scores[x]
    otu_id = otu[otu.find('_')+1:]
    d.update(c)
    q = insert_record.format(**d)
    db.execute(q, (otu_id,))
    
###
# Top level function: initialize the workspace directory, run the app

def classify(db, args):
    init_workspace(args)
    run_classifier(args)
    import_results(db, args)
    
###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    
    args = init_api(
        desc = "Use the RDP classifier to determine taxonomic categories for OTUs.",
        specs = [
            ('workspace',    { 'metavar': 'dir', 'help' : 'working directory', 'default' : 'taxonomy' } ),
            ('directory',    { 'metavar': 'dir', 'help' : 'directory containing mapped OTUs', 'default' : 'map' } ),
        ]
    )
        
    db = sqlite3.connect(args.dbname)
    record_metadata(db, 'start', ' '.join(sys.argv[1:]))

    # the taxonomy table has two columns for each level -- one is a foreign key for the
    # table of names for that level, the other is the probability of the classification

    try:
        taxonomy_spec = list(map(lambda x: (x+'_id', 'foreign', x), levels))
        taxonomy_spec += list(map(lambda x: ('p_'+x, 'REAL'), levels))
        init_table(db, 'taxonomy', 'taxonomy_id', taxonomy_spec, args.force)
    except Exception as err:
        print('Error while initializing output tables:', err)
        argparse.ArgumentParser.exit(1, 'Script aborted')

    classify(db, args)
    record_metadata(db, 'end', '')

    db.commit()
