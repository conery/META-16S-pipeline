#! /usr/bin/env python3

# Determine the taxonomic classification of a set of sequences

# John Conery / Kevin Xu Junjie
# University of Oregon
# 2014-06-05

import sqlite3
import argparse
import os
import os.path
from string import punctuation
import sys

# These global vars define filenames and other items shared by two or more steps

table_name = 'taxonomy'
input_file = 'otus.fasta'
output_file = 'classifications.txt'

# Path to the RDP classifier

# classifier_path = '~/Applications/rdp_classifier_2.7/dist/classifier.jar'
classifier_path = os.path.join(sys.path[0], 'rdp_classifier_2.7/dist/classifier.jar')

# This list is used to construct queries and dictionary keys used by the results parser.
# Notes
# **  'order' has an 'x' on the end -- it's a keyword in SQL and can't be used as a table name :-(
# ** RDP (v2.7) does not assign at the species level; the -f fixrank option generates these 6 ranks

levels = ['domain', 'phylum', 'class', 'orderx', 'family', 'genus']

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
# Initialize (or re-initialize) the tables that will hold the output from the classifier.

find_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{}"'.format(table_name)
drop_table = 'DROP TABLE {}'.format(table_name)
col_specs = ', '.join(map(lambda x: x + '_id INTEGER REFERENCES ' + x, levels))
col_specs += ', ' + ', '.join(map(lambda x: 'p_' + x + ' REAL', levels))
create_table = 'CREATE TABLE {} ( {}_id INTEGER PRIMARY KEY AUTOINCREMENT, otu_id INTEGER NOT NULL REFERENCES otu, {} )'.format(table_name, table_name, col_specs)
    
def prepare_tables(db, args):
    """
    Return True if the results table is initialized and ready to accept values.  If the
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
# Run the app

def run_classifier(args):
    cmnd = 'java -jar ' + classifier_path
    cmnd += ' classify ' + os.path.join(args.directory, input_file)
    cmnd += ' -f fixrank'
    cmnd += ' -o ' + os.path.join(args.workspace, output_file)
    print(cmnd)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'exec', cmnd))
    db.commit()
    res = os.system(cmnd)

###
# Create a map that associates barcodes with barcode IDs

# fetch_barcode_ids = 'SELECT barcode_id, sequence FROM barcodes'
#
# def barcode_map(db):
#     bc = { }
#     for bcid, sequence in db.execute(fetch_barcode_ids):
#         bc[sequence] = bcid
#     return bc

###
# Populate the table.  The function uses 7 maps, one for each taxonomic level.  As new
# names are discovered they are added to the map with a unique integer key.  After
# reading all the assignments the maps are written as new tables in the database.

def import_results(db, args):
    # inits = [s[0] for s in levels]
    # names = dict(zip(inits, [dict() for i in range(len(inits))]))
    names = dict(zip(levels, [dict() for i in range(len(levels))]))
    for line in open(os.path.join(args.workspace, output_file)):
        parts = dict(zip(levels, [None]*len(levels)))
        scores = dict(zip(levels, [0]*len(levels)))
        recs = line.split('\t')
        sanitize(recs)
        for i in range(2, len(recs), 3):
            parts[recs[i+1]] = recs[i]
            scores[recs[i+1]] = recs[i+2]
        update_names(names, parts)
        save_record(recs[0], names, parts, scores)
    save_tax_tables(names)
    
###
# Helper function for import_results -- create and fill the tables for each level
# of the taxonomy

def save_tax_tables(names):
    # tmap = dict(zip([s[0] for s in levels], levels))        # map initial to category
    for x in levels:
        db.execute('DROP TABLE IF EXISTS {}'.format(x))
        db.execute('CREATE TABLE {} ( {}_id INTEGER PRIMARY KEY, name TEXT)'.format(x,x))
    for x, values in names.items():
        for name, item_id in values.items():
            # q = 'INSERT INTO {} ({}_id, name) VALUES (?,?)'.format(tmap[x],tmap[x])
            q = 'INSERT INTO {} ({}_id, name) VALUES (?,?)'.format(x,x)
            db.execute(q, (item_id, name))
    
###
# Helper function for import_results -- iterate over the parts of a taxonomy string
# to update the name tables

def update_names(names, parts):
    for p, x in parts.items():          # p is a level (domain, phylum, etc), x is the name at that level
        if x is None or len(x) == 0: continue
        d = names[p]
        if x not in d:
            d[x] = len(d)+1

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
# Helper function -- change 'order' to 'orderx' (order is a keyword in SQL), strip
# quotes from names, etc

def sanitize(recs):
    for i in range(len(recs)):
        if recs[i] == 'order':
            recs[i] = 'orderx'
        recs[i] = recs[i].strip(punctuation)

###
# Top level function: initialize the workspace directory, run the app

def classify(db, args):
    init_workspace(args)
    run_classifier(args)
    import_results(db, args)

###
# Parse the command line arguments, call the top level function...

def init_api():
    parser = argparse.ArgumentParser(
        description="Use the RDP classifier to determine taxonomic categories for OTUs",
    )
    parser.add_argument('dbname', help='the name of the SQLite database file')
    parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
    parser.add_argument('-w', '--workspace', help='working directory', default='taxonomy')
    parser.add_argument('-d', '--directory', help='directory with input file', default='map')
    return parser.parse_args()
    
###
# Parse the command line arguments, call the top level function...
    
if __name__ == "__main__":
    args = init_api()
    
    db = sqlite3.connect(args.dbname)
        
    if not prepare_tables(db, args):
        argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')

    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
    classify(db, args)
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )

    db.commit()
