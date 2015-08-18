# FASTA sequence class, comparable to the FASTQ class (also part of this project).  

class FASTA:
	"""
	Simple FASTA sequence class
	"""
	def __init__(self, defline, sequence):
		self._def = defline
		self._seq = sequence
	
	def defline(self):
		"Return the defline string for this object."
		return self._def
		
	def sequence(self):
		"Return the sequence string for this object."
		return self._seq
			
	def __repr__(self):
		"""
		Return a 2-line string suitable for printing in a FASTA file.
		"""
		return '\n'.join(['>' + self._def, self._seq])

class FASTAReader:
	"""
	Create a FASTAReader object to read sequences from a file.
	"""
	def __init__(self, filename):
		self._f = open(filename)
		self._buf = self._f.readline()
		# shouldn't happen, but.... skip any blank lines at the front of the file
		while len(self._buf) == 1:
			self._buf = self._f.readline()

	def read(self):
		"""Return the next sequence from this file, or None if there are no more sequences"""
		if len(self._buf) == 0:									# _buf has defline or empty line if at the end of the file
			return None
		if self._buf[0] == '>':
			defline = self._buf.strip('>; \n')					# don't save the > at front or semicolon at end
		else:
			raise Exception("Defline does not start with '>'")
		sequence = ''
		self._buf = self._f.readline()
		while len(self._buf) > 0 and self._buf[0] != '>':		# when loop exits _buf has the first line of the next sequence
			sequence += self._buf.strip()
			self._buf = self._f.readline()
		return FASTA(defline, sequence)
