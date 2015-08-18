#! /usr/bin/env Rscript --vanilla

# Compute divsersity statistics for a set of 16S rRNA samples

# John Conery
# University of Oregon
# 2014-10-16

library(RSQLite)
library(vegan)

# Parse command line arguments -- get database name from first position, check
# to see if any options are specified.  Return a list of named options.

initialize <- function(args)
{
	filename <- args[1]
	if (is.na(filename))
	{
		message('Usage: diversity.R dbname [options]')
		quit(status = 1)
	}
	
	list(dbname=filename, force=is.element('--force', args))
}

# Make sure the file exists, open the database.  If the --force option is 
# specified remove any old table.

connect <- function(args)
{
	fn <- args$dbname
	
	if (! is.element(fn, Sys.glob('*.db')))
	{
		message('Database not found: ', fn)
		quit(status = 1)
	}

	db <- dbConnect(SQLite(), dbname=fn)
	
	if (args$force)
		for (tbl in c('alpha', 'evenness', 'beta'))
			if (dbExistsTable(db, tbl))
				dbRemoveTable(db, tbl)				
	
	db
}

# Top level starts here:  parse the command line arguments, open the database

args <- initialize(commandArgs(TRUE))
db <- connect(args)

# Read the abundance table (which has the main data) and samples table (which has the names of the experiments)

abundance <- dbReadTable(db, "abundance")
samples <- dbReadTable(db, "samples")

# Since the data is sorted by experiment, then by group, this command makes an n x m matrix with
# one row for each taxonomic group and one column for each sample:

m = matrix(abundance$n, nrow(abundance)/nrow(samples), nrow(samples))

# Assign column names using sample names:

colnames(m) <- samples$name

# When computing diversity statistics pass the transpose of the matrix since Vegan expects
# a frame with one row per sample and one column per taxononmic group

message('α diversity')
H <- diversity(t(m))
dbWriteTable(db, "alpha", data.frame(H))

message('evenness')
J <- H / log(specnumber(t(m)))
dbWriteTable(db, "evenness", data.frame(J))

message('β diversity')
beta <- vegdist(t(m), binary=TRUE)
dbWriteTable(db, "beta", data.frame(beta = matrix(beta)))


