#!/usr/bin/env python


# naive substring finding of subsequence p in larger sequence t
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences
    

# reverse complement input string
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t
    
    
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
    
    
# naive string comparison of subsequence p and its reverse complement
# with a larger sequence t
def naive_with_rc(p, t):
    occurrences = naive(p, t)
    p_rc = reverseComplement(p)
    if p_rc != p:
        occurrences += naive(p_rc, t)
        
    return occurrences
    

# naive substring finding of subsequence p in larger sequence t,
# allowing for up to 2 mismatches
def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        numMismatch = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                numMismatch += 1
                if numMismatch > 2:
                    break
        if numMismatch <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences
    
    
# section 1
lambdaGenome = readGenome('lambda_virus.fa')

AGGT_occur = naive_with_rc('AGGT', lambdaGenome)
print "num AGGT in lambda", len(AGGT_occur)

TTAA_occur = naive_with_rc('TTAA', lambdaGenome)
print "num TTAA in lambda:", len(TTAA_occur)

ACTAAGT_occur = naive_with_rc('ACTAAGT', lambdaGenome)
print "leftmost ACTAAGT in lambda:", min(ACTAAGT_occur)

AGTCGA_occur = naive_with_rc('AGTCGA', lambdaGenome)
print "leftmost AGTCGA in lambda:", min(AGTCGA_occur)

TTCAAGCC_occur = naive_2mm('TTCAAGCC', lambdaGenome)
print "num TTCAAGCC in lambda with <=2 mismatch:", len(TTCAAGCC_occur)

AGGAGGTT_occur = naive_2mm('AGGAGGTT', lambdaGenome)
print "leftmost AGGAGGTT in lambda with <=2 mismatch:", min(AGGAGGTT_occur)

# section 2

# returns an array containing the mean quality score at each cycle number
def qualByPosition(quals, readLength=100):
    qualScores = [0] * readLength
    totals = [0] * readLength
    
    for qual in quals:
        for i in range(len(qual)):
            qualScores[i] += ord(qual[i])-33
            totals[i] += 1
            
    for i in range(len(qualScores)):
        if totals[i] > 0:
            qualScores[i] /= float(totals[i])
        
    return qualScores
    
seqs, quals = readFastq('ERR037900_1.first1000.fastq')
qualScoresByPos = qualByPosition(quals)
minQual = min(qualScoresByPos)
print qualScoresByPos, minQual, qualScoresByPos.index(minQual)