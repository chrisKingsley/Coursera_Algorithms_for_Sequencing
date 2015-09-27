#!/usr/bin/env python

import sys

# read genome sequence from fasta file
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

    
# read seqs and quality scores from fastq file
def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities
    
    
def editDistanceSubstring(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the minimum value in the last row of the matrix
    return min(D[-1])
    

# chr1_seq = readGenome('../Homework2/chr1.GRCh38.excerpt.fasta')
# print editDistanceSubstring('GCTGATCGATCGTACG', chr1_seq)
# print editDistanceSubstring('GATTTACCAGATTGAG', chr1_seq)


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match
    
    

def printAllOverlaps(seqs, min_length):
    kmerHash = dict()
    numRead1WithOverlap, numPairs = 0, 0
    
    # populate dict where kmers of min_length point to a set of
    # reads which contain them
    for seq in seqs:
        for start in range(len(seq) - min_length + 1):
            kmer = seq[start:(start+min_length)]
            if kmer in kmerHash:
                kmerHash.get(kmer).add(seq)
            else:
                kmerHash[kmer] = { seq }
    
    # find the number of pairs with overlaps of at least min_length
    # and the number of reads that have one or more suffix overlaps
    for read1 in seqs:
        suffix = read1[ len(read1)-min_length: ]
        read1Overlap = False
        
        for read2 in kmerHash[ suffix ]:
            if read1!=read2:
                suffixLength = overlap(read1, read2, min_length)
                if suffixLength >= min_length:
                    numPairs += 1
                    read1Overlap = True
        if read1Overlap:
            numRead1WithOverlap += 1
                    
    print numPairs, numRead1WithOverlap      
    

min_length = 30
seqs, quals = readFastq('ERR266411_1.for_asm.fastq')
printAllOverlaps(seqs, min_length)
