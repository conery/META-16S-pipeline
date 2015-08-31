#! /usr/bin/env Rscript --vanilla

# Produce a stacked bar chart for set of 16S rRNA samples

# John Conery
# University of Oregon
# 2014-10-23

library(RSQLite)
library(ggplot2)

# Parse command line arguments -- get database name from first position, check
# to see if any options are specified.  Return a list of named options.

initialize <- function(args)
{
	dbfile <- args[1]
	outfile <- args[2]

	if (is.na(dbfile) | is.na(outfile))
	{
		message('Usage: barchart.R dbname X.pdf')
		quit(status = 1)
	}
		
	list(dbname=dbfile, filename=outfile)
}

# Make sure the file exists, open the database.  Make sure the database has a 
# table named abundance and that the group ID in that table is 'phylum_id'

connect <- function(args)
{
	fn <- args$dbname
	
	if (! is.element(fn, Sys.glob('*.db')))
	{
		message('Database not found: ', fn)
		quit(status = 1)
	}

	db <- dbConnect(SQLite(), dbname=fn)
	
	if (! (dbExistsTable(db, 'abundance') & is.element('phylum_id', dbListFields(db, 'abundance'))) )
	{
		message('Database does not have an abundance table grouped by phylum')
		quit(status = 1)
	}
	
	db
}

# Top level starts here:  parse the command line arguments, open the database

args <- initialize(commandArgs(TRUE))
db <- connect(args)

# Create a "long format" data frame.  Each row will have a sample name, a phylum
# name, and the count of the number of individuals.

f <- dbGetQuery(db, "select samples.name as sample, phylum.name as phylum, n from samples join abundance using (sample_id) join phylum using (phylum_id)")

# Make the bar chart

ggplot(f, aes(x = sample, y = n, fill = phylum)) + 
  geom_bar(stat="identity", colour="black") + 
  guides(fill=guide_legend(reverse=FALSE))
  
ggsave(file=args$filename)
