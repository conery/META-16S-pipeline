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
#   desc:         script description to print as part of the help message
#   specs:        a list of tuples with argument names and argument specs
#   epi:          a second string to print in the help message
#   with_limits:  early stages allow limits on the sample names or number of sequences

import argparse

def init_api(desc, specs, epi=None, with_limits=False):
    parser = argparse.ArgumentParser(description=desc, epilog=epi)
    parser.add_argument('dbname', help='the name of the project database file')
    parser.add_argument('--force', action='store_true', help='[GUI:flag::T] replace existing data')
    if with_limits:
        parser.add_argument('--noimport', action='store_true', help="test only, don't import results")
        parser.add_argument('--limit', metavar='N', type=int, help="import only N results per sample")
        parser.add_argument('--sample', metavar='id', help="use data from this sample only")
    for x in specs:
        parser.add_argument('--'+x[0], **x[1])
        
    args = parser.parse_args()
    return args

# Scripts early in the pipeline can restrict operations to a specified sample
# (e.g. when testing the pipeline). This function returns a list of records from
# the 'samples' table.  If a sample name is specfied on the command line
# use it, otherwise fetch all sample names.
#
# The return value is a list of records with sample ID, sample name, and FASTQ 
# file names for the specified samples (sample ID has to be the first column).

def sample_list(db, args):
    "Return a list of samples to process"
    sql = 'SELECT sample_id, name, r1_file, r2_file FROM samples'
    if isinstance(args.sample, str):
        sql += ' WHERE name = "{}"'.format(args.sample)
    return db.execute(sql).fetchall()

###
# Initialize the directory where intermediate work products will be stored.

import os

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
# Initialize the table that holds results of a stage.  To make sure users don't
# accidentally lose existing data this function makes sure the script won't
# overwrite current records unless the --force option is specified on the
# command line.  
#
# If there are any problems an exception is raised.  If the function returns
# the script can assume the table was initialized.

# Cases:
# (1) the output table doesn't exist; create it
# (2) data from this run won't conflict with existing data (the command line)
#      includes --sample X and no data for X is in the table)
# (3) there is conflicting data, --force is not specified; raise exception
# (4) remove conflicting data from the table

# Parameters:

def init_table(db, name, key, cols, replace, sample_name=None):
    """Initialize the output table for this stage in the pipeline.  Parameters:
           db           reference to database
           name         output table name
           key          name of primary key column
           cols         an array of tuples with specs for remaining columns
           replace      flag that specifies whether --force is on the command line
           sample_name  (optional) sample to be processed
    """
    
    def fetch_table_spec(db, name):
        sql = 'SELECT sql FROM sqlite_master WHERE name = "{x}"'.format(x=name)
        res = db.execute(sql).fetchall()
        return res[0][0] if res else False
    
    def col_spec(tup):
        'Convert a tuple into a SQL column spec'
        if tup[1] == 'foreign':
            return '%s INTEGER NOT NULL REFERENCES %s' % (tup[0], tup[2])
        else:
            return '%s %s' % (tup[0], tup[1])
    
    def create_table_query():
        'Combine the column specs, return the CREATE TABLE query'
        template = 'CREATE TABLE {tbl} ({key_col} INTEGER PRIMARY KEY AUTOINCREMENT, {cols})'
        s = ', '.join(map(col_spec, cols))
        return template.format(tbl=name, key_col=key, cols=s)

    def constrained(template, sample_id):
        'Add sample ID constraints if the sample name was specified'
        sql = (template + ' FROM {tbl}').format(tbl=name)
        if sample_id is not None:
            sql += ' WHERE sample_id = {x}'.format(x=sample_id)
        return sql
    
    def count(sample_id):
        'Look in the current table, see if there are records for this sample'
        sql = constrained('SELECT count(*)', sample_id)
        return db.execute(sql).fetchall()[0][0]
    
    def clear(sample_id):
        'Delete previous records for this sample'
        sql = constrained('DELETE', sample_id)
        record_metadata(db, 'query', sql)
        db.execute(sql)
        
    # See if the table exists; if not, make it and return
    existing = fetch_table_spec(db, name)
    if not existing:
        sql = create_table_query()
        record_metadata(db, 'query', sql)
        db.execute(sql)
        return
    
    # If the table has a 'sample_id' column and the command line has a --sample
    # option map a sample name into a sample ID; if the name is invalid abort the script
    if existing.find('sample_id') and sample_name is not None:
        res = db.execute('SELECT sample_id FROM samples WHERE name = ?', (sample_name, )).fetchall()
        if len(res) == 0:
            raise Exception('unknown sample name: ' + sample_name)
        sample_id = res[0][0]
    else:
        sample_id = None
    
    # The 'count' function returns the number of conflicting records
    if count(sample_id) == 0:
        return 
    
    # If --force specified remove the old records before returning.
    if replace:
        clear(sample_id)
        return
    
    raise Exception('Conflicting records, use --force to replace old data') 


###
# Add a message to the log table.  Arguments are a reference to the database,
# the type of event and a longer description.  The commit option updates the database,
# useful if a script is about to run a command that might crash (and thus abort
# the script itself without saving the DB).

import sys
import os

def record_metadata(db, event, message, commit=False):
    app = os.path.basename(sys.argv[0])
    db.execute("INSERT INTO log VALUES (DATETIME('NOW'), ?, ?, ?)", (app, event, message))
    if commit:
        db.commit()


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

