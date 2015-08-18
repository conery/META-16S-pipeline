#! /usr/bin/env python3

# Print the abundance table in spreadsheet form, labeling rows by taxonomic units
# and columns by experiment names

# John Conery
# University of Oregon
# 2014-10-17

import sqlite3
import argparse
import sys

# Taxonomic groups, listed in order

groups = [ 'domain', 'phylum', 'class', 'order', 'family', 'genus' ]

# Map taxonomic groups to table names. Note that 'order' has an 'x' on the end -- 
# 'order' is a keyword in SQL and can't be used as a table name :-(

tax_tbl_name = { x : x for x in groups }
tax_tbl_name['order'] += 'x'
	
###
# Get the table description, figure out how it is grouped.  Stupid approach, just
# returns the first taxonomic table name found in the abundance table description.

fetch_desc = 'SELECT sql FROM sqlite_master WHERE name = "abundance"'

def find_group(db):
	schema = db.execute(fetch_desc).fetchall()[0]
	if schema:
		for tbl in tax_tbl_name:
			if tbl in schema[0]:
				return tbl
	return None

### 
# Get the list of experiment names.  Important: the names should be sorted by the
# same criteria used to fetch abundances.

def fetch_sample_names(db):
	return db.execute('SELECT name FROM samples ORDER BY name').fetchall()
	
###
# Create and exdecute a query that fetches the abundance table along with taxonomic category names

def fetch_abundance_table(db, grp):
	query = 'SELECT samples.name, n'
	for x in groups:
		query += ', {g}.name AS {g}'.format(g=tax_tbl_name[x])
		if x == grp:
			break
	query += ' FROM samples JOIN abundance USING (sample_id) JOIN taxonomy USING ({g}_id)'.format(g=tax_tbl_name[grp])
	for x in groups:
		query += ' JOIN {g} USING ({g}_id)'.format(g=tax_tbl_name[x])
		if x == grp:
			break
	query += ' GROUP BY sample_id, {g}_id'.format(g=tax_tbl_name[x])
	query += ' ORDER BY ' + tax_tbl_name[groups[0]]
	if grp != groups[0]:
		for x in groups[1:]:
			query += ', ' + tax_tbl_name[x]
			if x == grp:
				break
	return db.execute(query).fetchall()

###
# The first line should have blank columns for the taxonomic groups, then
# the names of the experiments

def print_header(nlevels, exp):
	pre = [''] * nlevels
	post = list(map(lambda x: x[0], exp))
	print(','.join(pre+post))
	
###
# Compare the previous taxonomy columns with the current record, return a
# new sequence with the groups that changed.

def tax_cols(rec, prev):
	a = []
	for i in range(nlevels):
		if rec[i+2] == prev[i]:
			a.append('')
		else:
			a.append(rec[i+2])
			prev[i] = rec[i+2]
	if rec[1] is None:
		a.append('')
	else:
		a.append(str(rec[1]))
	return a

def print_table(nlevels, ncols, recs):
	nrows = len(recs) // ncols
	prev = [None] * nlevels
	for i in range(nrows):
		# compare the groups in the current row with the previous set of groups
		# to create the first part of the line
		cur = recs[i*ncols]
		a = []
		for j in range(nlevels):
			if cur[j+2] == prev[j]:
				a.append('')
			else:
				a.append(cur[j+2])
				prev[j] = cur[j+2]
		# append the counts for the current row
		for j in range(ncols):
			cur = recs[i*ncols + j]
			a.append(str(cur[1]))
		# print the row
		print(','.join(a))

###
# Parse the command line arguments, call the top level function...

def init_api():
	parser = argparse.ArgumentParser(
		description="Print the abundance table in spreadsheet format",
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	return parser.parse_args()
	
###
# Parse the command line arguments, call the top level function...
	
if __name__ == "__main__":
	args = init_api()
	db = sqlite3.connect(args.dbname)
	grp = find_group(db)
	if grp is None:
		print(sys.argv[0], 'missing abundance table or unknown group')
		sys.exit(1)
		
	res = fetch_abundance_table(db, grp)
	samples = fetch_sample_names(db)
	nlevels = groups.index(grp) + 1
	
	print_header(nlevels, samples)
	print_table(nlevels, len(samples), res)
	
	# prev = res[0]
	# for x in res:
	# 	print(','.join())

## First version, saves taxonomic group names -- save this for script that prints data in human
## readable forms

# def create_abundance_table(db, args):
# 	query = 'CREATE TABLE {t} AS SELECT barcode_id, '.format(t=table_name)
# 	for x in groups:
# 		query += '{g}.name AS {g}, '.format(g=tax_tbl_name[x])
# 		if x == args.group:
# 			break
# 	query += 'count(*) as otus, sum(count) AS n'
# 	query += ' FROM otu JOIN taxonomy USING (otu_id)'
# 	for x in groups:
# 		query += ' JOIN {g} USING ({g}_id)'.format(g=tax_tbl_name[x])
# 		if x == args.group:
# 			break
# 	query += ' GROUP BY barcode_id, {g}.name'.format(g=tax_tbl_name[args.group])
# 	query += ' ORDER BY barcode_id, ' + tax_tbl_name[groups[0]]
# 	if args.group != groups[0]:
# 		for x in groups[1:]:
# 			query += ', ' + tax_tbl_name[x]
# 			if x == args.group:
# 				break
# 	print(query)
# 	db.execute(query)


## Start on a version that uses outer joins -- has glitches in connecting to taxonomy table...

# def create_abundance_table(db, args):
# 	query = 'CREATE TABLE {t} AS SELECT experiment, {g}.name AS {g}'.format(t=table_name, g=tax_tbl_name[args.group])
# 	for x in groups:
# 		if x == args.group:
# 			break
# 		query += ', {g}.name AS {g}'.format(g=tax_tbl_name[x])
# 	query += ', count(*) as otus, sum(count) AS n'
# 	query += ' FROM barcodes JOIN {g}'.format(g=tax_tbl_name[args.group])
# 	query += ' LEFT OUTER JOIN otu USING (barcode_id, {g}_id) JOIN taxonomy USING (otu_id)'.format(g=tax_tbl_name[args.group])
# 	for x in groups:
# 		if x == args.group:
# 			break
# 		query += ' JOIN {g} USING ({g}_id)'.format(g=tax_tbl_name[x])
# 	query += ' GROUP BY barcode_id, {g}.name'.format(g=tax_tbl_name[args.group])
# 	query += ' ORDER BY barcode_id, ' + tax_tbl_name[groups[0]]
# 	if args.group != groups[0]:
# 		for x in groups[1:]:
# 			query += ', ' + tax_tbl_name[x]
# 			if x == args.group:
# 				break
# 	print(query)
# 	db.execute(query)
