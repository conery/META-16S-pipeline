#! /usr/bin/env python3

#
# demonstrate how to use the FASTQ class
#

from FASTQ import *

def print_info(seq):
    print()
    print('defline:  ', seq.defline())
    print('sequence: ', seq.sequence())
    print('quality:  ', seq.quality())
    print('blob:     ', seq.blob())

# Create a FASTQ object from a string

seq1 = FASTQ('@NS500451\nNCAAGNGCCA\n+\n#AAAA#FFFF')

# The _repr_ method returns a string suitable for printing in a FASTQ file:

print('Test Sequence 1:')
print(seq1)

# Print the attributes:

print_info(seq1)

# Compress the sequence, print the attributes again:

seq1.pack()
print_info(seq1)

# We should be able to restore the original attributes:

seq1.unpack()
print_info(seq1)

# Create another FASTQ object, this time from a byte array (e.g. a blob
# written into a database).  The sequence will have an empty defline, so
# we need to add one:

seq2 = FASTQ(b'\xc8P\x1e\x1e\xb4\xc8\xb9UU#')
seq2._def = 'restored sequence'
seq2.unpack()

print('\nTest Sequence 2:')
print(seq2)
print()

# This test shows how to read sequences from a file using the FASTQReader class.
# Since these are real sequences with actual deflines we can also demo the
# defline parser classes.

fn = 'test.fastq'
print('Sequences in file', fn)

with FASTQReader(fn) as file:
    for seq in file:
        print()
        print(seq.defline())
        print(len(seq), "bp")
        seq.parse_defline()
        print('lane:tile:x:y', seq.unique_defline())
        print('filtered: ', seq.filtered())

# The readseq method reads one read at a time.  If readline or another method
# gets a partial read a call to readseq will skip ahead to the start of the 
# next read.

file = FASTQReader(fn)
seq1 = file.readseq()
def2 = file.readline().strip()
seq3 = file.readseq()

print('\nThree sequence deflines:')
print(seq1.defline())
print(def2)
print(seq3.defline())

# The next call should return None

print(file.readseq())



