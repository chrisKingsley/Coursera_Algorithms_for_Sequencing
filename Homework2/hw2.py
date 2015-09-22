#!/usr/bin/env python

from bm_preproc import *
from kmer_index import *

# read genome sequence from fasta file
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    
    
# naive substring finding of subsequence p in larger sequence t
def naive_with_counts(p, t):
    occurrences = []
    numAligns, numCharCompares = 0, 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        numAligns += 1
        for j in range(len(p)):  # loop over characters
            numCharCompares += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return (occurrences, numAligns, numCharCompares)
    
    
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences
    
    
def boyer_moore_with_counts(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    numAligns, numCharCompares = 0, 0
    while i < len(t) - len(p) + 1:
        shift = 1
        numAligns += 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            numCharCompares += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return (occurrences, numAligns, numCharCompares)

ch1_portion = readGenome('chr1.GRCh38.excerpt.fasta');
aluSeq = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'

# print "Naive matching:",
# print naive_with_counts(aluSeq, ch1_portion)

# p_bm = BoyerMoore(aluSeq)
# print "Boyer Moore:",
# print boyer_moore_with_counts(aluSeq, p_bm, ch1_portion)


def approx_match(p, t, n):
    segmentLength = int(round(len(p) / (n+1)))
    allMatches = set()
    numIndexHits = 0
    for i in range(n+1):
        start = i*segmentLength
        end = min((i+1)*segmentLength, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        numIndexHits += len(matches)
        
        for m in matches:
            if m<start or m-start+len(p)>len(t):
                continue
                
            mismatches = 0
            for j in range(0, start):
                if not p[j]==t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j]==t[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
                        
            if mismatches <= n:
                allMatches.add(m-start)
                
    return list(allMatches), numIndexHits

    
short_alu = 'GGCGCGGTGGCTCACGCCTGTAAT'
# matches, numIndexHits = approx_match(short_alu, ch1_portion, 2)
# print matches, len(matches), numIndexHits



class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits
        
        
def approxMatchSubSeq(p, t, subIdx, n):
    segmentLength = int(round(len(p) / (n+1)))
    print segmentLength
    allMatches = set()
    numIndexHits = 0
    for i in range(n+1):
        matches = subIdx.query(p[i:])
        print p[i:], matches
        numIndexHits += len(matches)
        
    return numIndexHits
    
    
subIdx = SubseqIndex(ch1_portion, 8, 3)
print approxMatchSubSeq(short_alu, ch1_portion, subIdx, 2)

