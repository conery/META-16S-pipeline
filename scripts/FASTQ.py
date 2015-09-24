# FASTQ sequence class

# John Conery / Kevin Xu
# University of Oregon
# craeted Dec 2014, updated Aug 2015

# An object of the FASTQ class represents a single sequence read from a FASTQ
# format sequence file.  Objects may be compressed or uncompressed.  To make an
# uncompressed sequence pass the constructor a single string with 4 lines of text 
# read from a file; the constructor will split the lines and save the 1st line 
# as a defline attribute, the 2nd as the sequence, and the 4th as the quality string.  

# To compress a sequence call the pack method to compress both the sequence and 
# quality string into a single byte array.  The sequence and quality attributes 
# will be replaced by the byte array, which is called 'blob'.  To uncompress a 
# packed object call unpack, which restores the sequence and quality attributes.  

# A new object can also be created by passing a blob (e.g. a byte string read 
# from a database) to the constructor.  When created this way the defline is
# empty; it can be assigned later by assigning to self._def.

# Use 'getter' methods to access components of a sequence:
#   defline     the defline if the sequence was read from a file, else ''
#   sequence    the sequence letters if the sequence is uncompressed, else None
#   quality     the quality letters if the sequence is uncompressed, else None
#   blob        the compressed byte string, or None if the squence is uncompressed
#   filtered    True if the defline has a 'Y' in the filtered field

# It's also possible to access the parts of a sequence directly:
#   _def      the original defline, as read from the file
#   _seq      the sequence letters
#   _qual     the quality characters
#   _blob     a byte string for a compressed sequence
# Note that _blob is not defined until a sequence is compressed, and when it is
# compressed the _seq and _qual attributes are deleted.

# A method named parse_defline will define additional attributes extracted from
# the defline (assuming it's from an Illumina file): _instrument, _run_id, 
# _flowcell_id, _lane, _tile, _x, _y, _pair_end, _filtered, _control_bits,
# and _index.  

# Call a method named unique_defline to get a string created from lane, tile, and 
# (x,y) coords, which will be sufficient to distinguish this read from all others
# in its data set.  

# The unique_id method is the same as unique_defline but it returns an integer.

import re

class FASTQ:
    "A DNA string with quality scores"
    
    def __init__(self, arg):
        "arg can be a string containing 4 lines read from a FASTQ file or a byte array containing a compressed sequence"
        if isinstance(arg, str):
            self._def, self._seq, dummy, self._qual = arg.split('\n')
        else:
            self._def = ''
            self._blob = arg
        self._parsed = False
    
    def defline(self):
        "Return the defline read from a FASTQ file."
        return self._def

    def sequence(self):
        "Return the sequence string for this object, or None if the object is compressed."
        return self._seq if hasattr(self, '_seq') else None
    
    def quality(self):
        "Return the quality string for this object, or None if the object is compressed."
        return self._qual if hasattr(self, '_qual') else None
    
    def blob(self):
        "Return the packed byte string if the object is compressed, otherwise None."
        return self._blob if hasattr(self, '_blob') else None
        
    def parse_defline(self):
        "split the Illumina format defline into components, save attributes for each component"
        fields = self._def.split()
        self._instrument, self._run_id, self._flowcell_id, self._lane, self._tile, self._x, self._y = fields[0].split(':')
        self._pair_end, self._filtered, self._control_bits, self._index = fields[1].split(':')
        self._parsed = True
        
    def set_defline(self, instrument, run_id, flowcell_id, lane, tile, x, y, pair_end, filtered, control, index):
        part1 = ':'.join(map(str, [instrument,run_id,flowcell_id,lane,tile,x,y]))
        part2 = ':'.join(map(str, [pair_end, filtered, control, index]))
        self._def = '@' + part1 + ' ' + part2
        
    def filtered(self):
        "Return True if the defline contains ':Y:'"
        return re.search(r':Y:', self._def) is not None
        
    def unique_defline(self):
        "Return a string with the unique identifying parts of this sequence"
        if self._parsed:
            return ':'.join([self._lane, self._tile, self._x, self._y])
        else:
            return None
    
    def unique_id(self):
        "Return a unique identifying integer for this sequence"
        raise Exception('TBD')

    int2char = 'ACTGN'
    char2int = {'A': 0, 'C':1, 'T':2, 'G':3, 'N':4}
    qbase = ord('#')
    
    def pack(self):
        "Use a lossless encoding scheme to compress the sequence and quality strings into a single byte array."
        self._blob = bytearray()
        for i in range(len(self._seq)):
            self._blob.append(50*FASTQ.char2int[self._seq[i]] + ord(self._qual[i])-FASTQ.qbase)
        self._blob = bytes(self._blob)
        delattr(self, '_seq')
        delattr(self, '_qual')
        
    def unpack(self):
        "Uncompress a compressed object."
        self._seq = ''
        self._qual = ''
        if self._blob is not None:
            for byte in self._blob:
                s, q = divmod(byte, 50)
                self._seq += FASTQ.int2char[s]
                self._qual += chr(FASTQ.qbase + q)
        delattr(self, '_blob')
    
    def __repr__(self):
        """
        The representation of an unpacked sequence is a 4-line string suitable for printing in 
        a FASTQ file.
        """
        a = [self._def]
        b = [self._seq, '+', self._qual] if hasattr(self, '_seq') else [str(self._blob)]
        return '\n'.join(a+b)
        
    def __len__(self):
        if self._seq:
            return len(self._seq)
        else:
            return len(self._blob)

    @staticmethod
    def parse_defline(s):
        """
        Remove redundant characters from an Illumina format defline, return a new defline
        formed from the lane, tile number, and x and y coordinates.
        """
        x = s.split(':')
        return ':'.join(x[3:7])
    
# A FASTQReader is a type of text file.  It adds a readseq method that returns a FASTQ object
# for the next sequence in the file.  It also implements the iterator pattern so users can
# iterate over the file.

# The readsseq method checks to make sure the sequence starts with a line that has a '@'
# in the first column.  If not, it skips ahead to the next line that starts with '@'.

import io

class FASTQReader(io.TextIOWrapper):
    def __init__(self, fn):
        "Make a new FASTQReader for sequences in file 'fn'"
        super().__init__(open(fn, 'rb'))

    def __iter__(self):
        return self
    
    def __next__(self):
        res = self.readseq()
        if res is None:
            raise StopIteration
        return res
        
    def readseq(self):
        "Return a FASTQ object for the next sequence in the file."
        res = self.readline()
        while len(res) > 0 and res[0] != '@':
            res = self.readline()
        if len(res) == 0:
            return None
        for i in range(3):
            res += self.readline()
        return FASTQ(res.strip())



