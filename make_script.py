#! /usr/bin/env python3

# Generate a script that will run the entire pipeline, based on log messages
# in a database.  Use the most recent invocation of each script.

# John Conery
# University of Oregon
# 2014-04-15

import sqlite3
import argparse

fetch_events = 'SELECT max(time), message FROM log WHERE event = "start" GROUP BY script ORDER BY time'

aciss_prolog = """
#!/bin/bash -l

#PBS -q generic
#PBS -m abe
#PBS -M conery@uoregon.edu

module load python/3.3.4
module load pandaseq
module load usearch

# cd /tmp
# mkdir 16S

"""

aciss_epilog = """
# cp stuff back to home dir
"""
	
def print_script(db, aciss):
	if aciss:
		print(aciss_prolog)
	else:
		print('#! /usr/bin/env bash')
	for time, script in db.execute(fetch_events).fetchall():
		print(script)
	if aciss:
		print(aciss_epilog)
	
### 
# Set up command line arguments

def init_api():
	parser = argparse.ArgumentParser(
		description="Use log messages in a PIP/NGS database to generate a script that will run the full pipeline.",
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('-a', '--aciss', action='store_true', help='include commands to run script on ACISS')
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()
	db = sqlite3.connect(args.dbname)
	print_script(db, args.aciss)
