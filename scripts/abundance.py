#! /usr/bin/env python3

# Summarize OTU classifications by a specified grouping

# John Conery
# University of Oregon
# 2014-10-03

import sqlite3
import argparse
import sys

# Output table name

table_name = 'abundance'

# Taxonomic groups, listed in order

groups = [ 'domain', 'phylum', 'class', 'order', 'family', 'genus' ]

# Map taxonomic groups to table names. Note that 'order' has an 'x' on the end -- 
# 'order' is a keyword in SQL and can't be used as a table name :-(

tax_tbl_name = { x : x for x in groups }
tax_tbl_name['order'] += 'x'

###
# Check to see if there is already an abundance table, delete it if --force specified
# on the command line

find_table = 'SELECT name FROM sqlite_master WHERE type = "table" AND name = "{}"'.format(table_name)
drop_table = 'DROP TABLE {}'.format(table_name)
	
def prepare_tables(db, args):
	"""
	Return True if the DB is ready to make a new table.  If the
	table exists already return False unless --force was specified on the command line.
	"""
	if db.execute(find_table).fetchall():
		if args.force:
			db.execute(drop_table)
		else:
			return False
	return True	 
	
###
# Top level function: create the abundance table by executing a query that
#  (a) does a cross join to create all combinations of the specified group and samples
#  (b) does a regular join to connect the taxonomy table to the OTU table
#  (c) does an outer join of the cross product with the OTUs and their specified group
# The query includes a 'group by' to collapse the rows according to the specified group

def create_abundance_table(db, args):
	query = 'CREATE TABLE {t} AS '.format(t=table_name)
	query += 'SELECT {g}.{g}_id AS {g}_id, samples.sample_id AS sample_id, sum(count) AS n '.format(g=tax_tbl_name[args.group])
	join1 = '{g} JOIN samples'.format(g=tax_tbl_name[args.group])
	join2 = 'otu JOIN taxonomy USING (otu_id)'
	query += 'FROM {} LEFT OUTER JOIN ({}) USING (sample_id, {}_id) '.format(join1, join2, tax_tbl_name[args.group])
	query += 'GROUP BY samples.sample_id, {g}.{g}_id'.format(g=tax_tbl_name[args.group])
	query += ' ORDER BY sample_id, {g}_id'.format(g=tax_tbl_name[args.group])
	print(query)
	db.execute(query)
	
	db.execute('UPDATE abundance SET n = 0 WHERE n IS NULL')

###
# Parse the command line arguments, call the top level function...

def init_api():
	parser = argparse.ArgumentParser(
		description="Create a table with the total number of OTUs for each combination of taxonomic group and experiment.",
		epilog="Possible arguments for --group:  " + ', '.join(groups)
	)
	parser.add_argument('dbname', help='the name of the SQLite database file')
	parser.add_argument('-f', '--force', action='store_true', help='re-initialize an existing table')
	parser.add_argument('-g', '--group', help='taxonomic level', default='genus')
	return parser.parse_args()
		
if __name__ == "__main__":
	args = init_api()

	db = sqlite3.connect(args.dbname)
		
	if not prepare_tables(db, args):
		argparse.ArgumentParser.exit(1, 'Table exists; use --force if you want to replace previous values')
		
	if args.group not in groups:
		print('Valid group names:', ', '.join(groups))
		argparse.ArgumentParser.exit(1, 'Unknown group: {}'.format(args.group))

	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'start', ' '.join(sys.argv)) )
	create_abundance_table(db, args)
	db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (sys.argv[0], 'end', '') )

	db.commit()

