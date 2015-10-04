#!/usr/bin/env python

import itertools, time


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
    
    
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match


def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = set()
    
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if len(shortest_sup)==0 or len(sup) < len(shortest_sup[0]):
            shortest_sup = [sup]  # found shorter superstring
        if len(shortest_sup)>0 and len(sup)==len(shortest_sup[0]):
            shortest_sup.append(sup)
            
    return set(shortest_sup)  # return shortest
    
# strings = [ 'CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT' ]
# shortest_sup = scs(strings) 
# print shortest_sup, len(list(shortest_sup)[0]), len(shortest_sup)


def pick_maximal_overlap(reads, k):
    reada, readb = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            reada, readb = a, b
            best_olen = olen
            
    return reada, readb, best_olen
        
     
def greedy_scs(reads, k):
    read_a, read_b, olen = pick_maximal_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
        
    return ''.join(reads)
    

def removeFromKmerDict(kmerDict, seq, min_length):
    for start in range(len(seq) - min_length + 1):
        kmer = seq[start:(start+min_length)]
        kmerSet = kmerDict.get(kmer)
        if seq in kmerSet:
            kmerSet.remove(seq)
        
def addToKmerDict(kmerDict, seq, min_length):
    for start in range(len(seq) - min_length + 1):
        kmer = seq[start:(start+min_length)]
        if kmer in kmerDict:
            kmerDict.get(kmer).add(seq)
        else:
            kmerDict[kmer] = { seq }
    
def getKmerDict(seqs, min_length):
    kmerDict = dict()
    numRead1WithOverlap, numPairs = 0, 0
    
    # populate dict where kmers of min_length point to a set of
    # reads which contain them
    for seq in seqs:
        addToKmerDict(kmerDict, seq, min_length)
                
    return kmerDict
    

    
def pick_maximal_overlap2(reads, kmerDict, k):
    reada, readb = None, None
    best_olen = 0
    for a in reads:
        suffix = a[ len(a)-k: ]
        for b in kmerDict[ suffix ]:
            if a!=b:
                olen = overlap(a, b, min_length=k)
                if olen > best_olen:
                    reada, readb = a, b
                    best_olen = olen
            
    return reada, readb, best_olen
    
def greedy_scs2(reads, kmerDict, k):
    t = time.time()
    read_a, read_b, olen = pick_maximal_overlap2(reads, kmerDict, k)
    
    while olen > 0:
        print len(reads), time.time() - t
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        
        removeFromKmerDict(kmerDict, read_a, k)
        removeFromKmerDict(kmerDict, read_b, k)
        addToKmerDict(kmerDict, read_a + read_b[olen:], k)
        
        t = time.time()
        read_a, read_b, olen = pick_maximal_overlap2(reads, kmerDict, k)
        
    return ''.join(reads)


k = 3
seqs, _ = readFastq('ads1_week4_reads.fq')
kmerDict = getKmerDict(seqs, k)
viralGenome = greedy_scs2(seqs, kmerDict, k)
print viralGenome, len(viralGenome)