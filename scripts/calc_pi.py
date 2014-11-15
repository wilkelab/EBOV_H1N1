#!/usr/bin/python
from Bio import SeqIO
import sys


def calc_distance(s1, s2):
    assert len(s1) == len(s2)
    d = 0.
    sites = 0
    allowed_chars = ('a', 'A', 'g', 'G', 'c', 'C', 't', 'T')
    for i in xrange(len(s1)):
        if s1[i] in allowed_chars and s2[i] in allowed_chars:
            sites += 1
            if s1[i] != s2[i]:
                d += 1.
    if sites < .9 * len(s1):
        print "Warning: removing >10% of data in comparison"
    return d/sites

def calc_pi(filename):
    handle = open(filename, "rU")
    seqs = []
    for record in SeqIO.parse(handle, "fasta") :
        seqs.append(record.seq)
    handle.close()

    d = 0.
    count = 0
    for i in xrange(len(seqs)):
        for j in xrange(i):
            d += calc_distance(seqs[i], seqs[j])
            count += 1
    print d/count

if len(sys.argv) != 2:
    print "Please provide alignment filename as command-line argument"
else:
    calc_pi(sys.argv[1])
