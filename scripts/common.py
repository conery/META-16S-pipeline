# A place to put code used by two or more scripts....

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

