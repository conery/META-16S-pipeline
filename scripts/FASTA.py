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
        return '\n'.join([self._def, self._seq])

    def __len__(self):
        return len(self._seq)

# A FASTAReader is a type of text file.  It adds a readseq method that returns a FASTA object
# for the next sequence in the file.  It also implements the iterator pattern so users can
# iterate over the file.

# The readsseq method checks to make sure the sequence starts with a line that has a '>'
# in the first column.  If not, it skips ahead to the next line that starts with '>'.  The
# method allows sequences to be spread across multiple lines.

import io
import sys

class FASTAReader(io.TextIOWrapper):
    def __init__(self, desc):
        "Make a new FASTAReader for sequences in file 'fn'"
        super().__init__(open(fn, 'rb'))
        # is there a way to read from stdin?  this code leads to an error ("underlying read
        # should have returned a bytes object") that was fixed by using 'rb' when opening
        # a file....
        # if isinstance(desc, str):
        #     super().__init__(open(desc, 'rb'))
        # else:
        #     super().__init__(sys.stdin)
        self._buffer = self.readline()             # initialize the buffer

    def __iter__(self):
        return self
    
    def __next__(self):
        res = self.readseq()
        if res is None:
            raise StopIteration
        return res

    def readseq(self):
        "Return a FASTA object for the next sequence in the file."
        while len(self._buffer) > 0 and self._buffer[0] != '>':
            self._buffer = self.readline()
        if len(self._buffer) == 0:
            return None
        defline = self._buffer.strip()
        sequence = ''
        self._buffer = self.readline()
        while len(self._buffer) > 0 and self._buffer[0] != '>':
            sequence += self._buffer.strip()
            self._buffer = self.readline()
        return FASTA(defline, sequence)
        