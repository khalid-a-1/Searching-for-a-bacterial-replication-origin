import sys
import random
import math


def getSequences(fname):
    """This function reads in the sequences from
    the given .fasta file."""
    seqs = []
    with open(fname, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seqs.append(line.rstrip())
    return seqs


def getCounts(motifs):
    """Return nt counts given a list of motif strings."""

    counts = []
    for i in range(len(motifs[0])):
        counts.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})

    for motif in motifs:
        
        motif = motif.upper()
        for i in range(len(motif)):
            counts[i][motif[i]] += 1

    return counts

def getConsensus(motifs):
    """Return the consensus sequence for a set of motifs."""

    counts = getCounts(motifs)

    consensus = ''
    for i in range(len(counts)):
        majority = 'A'
        for nt in counts[i]:
            if counts[i][nt] > counts[i][majority]:
                majority = nt
        consensus += majority

    return consensus

def getScore(motifs):
    """Return the score for a set of motifs."""

    counts = getCounts(motifs)
    consensus = getConsensus(motifs)

    t = len(motifs)
    score = 0
    for i in range(len(consensus)):
        nt = consensus[i]
        diff = t - counts[i][nt]
        score += diff

    return score

def probabilities(kmer, grid):
    """This function returns a probability matrix
    given a motif and it's profile."""
    pp = {'A': 0,'C': 1,'G': 2,'T': 3,}
    prob = 1
    for x,y in enumerate(kmer):
        prob *= grid[x][pp.get(y)]
    return prob

def probabableMotifs(sequence, profile):
    """This function return the most probable k-mer
    in sequence, given a profile."""
    k = len(profile)
    b_prob = 0
    b_kmer = None
    s =[]
    for b in range(len(sequence) - k+1):
        s.append(sequence[b:b+k])
    for i in s:
        p = probabilities(i, profile)
        if p > b_prob:
            b_kmer = i
            b_prob = p
    return b_kmer


def getProfile(motifs):
    """This function returns a profile from a set of motifs."""
    prof = []
    for i in zip(*motifs):
        l = {'A': 1.0, 'C': 1.0, 'G': 1.0, 'T': 1.0}
        for k in i:
            l[k] += 1
        for h in l:
            l[h] = l.get(h)/len(i)
        prof.append([l.get('A'), l.get('C'), l.get('G'), l.get('T')])
    return prof

def gibbs_sampler(dna, k, t, N):
    """This function implements randomized motif search
    using Gibbs Sampling and then returns the best set of
    motifs. """
    motifs = []
    for d in dna:
        ran = random.randint(0, k-1)
        motifs.append(d[ran:ran+k])
    best_motifs = motifs
    for i in range(N):
        i = random.randint(0, t-1)
        motifs.pop(i)
        profile = getProfile(motifs)
        motifs.insert(i, probabableMotifs(dna[i], profile))
        if getScore(motifs) < getScore(best_motifs):
            best_motifs = motifs
    return best_motifs


def getBestMotifSet(dna, k, t, N, numOfRepeats):
    """Returns the best motif set and best score after 
    a number of repeats of the motif search function"""

    bestMotifs = []
    bestScore = k * t
    
    for i in range(numOfRepeats):
        motifs = gibbs_sampler(dna, k, t, N)
        score = getScore(motifs)
        
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
    
    return bestMotifs, bestScore

def main():
    sequences = getSequences("upstream250.fasta")


    # d = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG','TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC','AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
    gibbs_output = gibbs_sampler(sequences, 20, 36, 1000)
    print(gibbs_output)
    
    best_motifs, bestScore = getBestMotifSet(sequences, 20, 36, 300, 2000)
    print(best_motifs, bestScore)

    consensus = getConsensus(best_motifs)
    print(consensus)

    motifs_from_tool = [
        'TTCGTGACCGACGTCCCCAG',
        'TTGGGGACTTCCGGCCCTAA',
        'GCCGGGACTTCAGGCCCTAT',
        'CATGGGACTTTCGGCCCTGT',
        'GAGGGGACTTTTGGCCACCG',
        'CCAGGGACCTAATTCCATAT',
        'TTGAGGACCTTCGGCCCCAC',
        'CTGGGGACCGAAGTCCCCGG',
        'TTAGGGACCATCGCCTCCTG',
        'TGGATGACTTACGGCCCTGA',
        'TTGGGGACTAAAGCCTCATG',
        'TCGGGGACTTCTGTCCCTAG',
        'TTGGGGACCATTGACCCTGT',
        'TTGAGGACCTAAGCCCGTTG',
        'CACGGGTCAAACGACCCTAG',
        'GGCGGGACGTAAGTCCCTAA',
        'GAAGTGACGAAAGACCCCAG',
        'CGGAGGACCTTTGGCCCTGC',
        'GTGGGGACCAACGCCCCTGG'
    ]

    print(getScore(motifs_from_tool))


main()
