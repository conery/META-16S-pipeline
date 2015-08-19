###
#
#  Shared utilities used by scripts in the META 16S Analysis Pipeline
#
###

# Each script requires the name of the project database as a positional argument
# and requires a --force argument to overwrite an existing table.  This function
# initializes an argument parser with common arguments, adds script-specific
# arguments, and calls the argument parser.
#
#   desc:   script description to print as part of the help message
#   specs:  a list of tuples with argument names and argument specs
#   epi:    a second string to print in the help message

import argparse

def init_api(desc, specs, epi=None):
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument('dbname', help='the name of the project database file')
    parser.add_argument('--force', action='store_true', help='[GUI:flag::T] re-initialize an existing database')
    
    for x in specs:
        parser.add_argument('--'+x[0], **x[1])
        
    args = parser.parse_args()
    return args
    
###
# Add a message to the log table.  Arguments are a reference to the database,
# the type of event and a longer description

import sys
import os

def record_metadata(db, event, message):
    app = os.path.basename(sys.argv[0])
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (app, event, message))
    

###
# Create a map that associates a defline with the sequence ID in the 
# panda table

fetch_from_panda = 'SELECT panda_id, defline FROM panda'

def defline_map(db, sample_id = None):
    dm = { }
    query = fetch_from_panda
    if sample_id:
        query += ' WHERE sample_id = {}'.format(sample_id)
    for pid, defline in db.execute(query):
        dm[defline] = pid
    return dm

