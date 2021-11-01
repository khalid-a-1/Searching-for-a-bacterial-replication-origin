import matplotlib.pyplot as plt

def readFile(fname):
    """
    Reads in .fasta file and returns the DNA sequence string.
    PARAMETERS:
        fname - file name
    RETURN:
        seq - DNA sequence string
    """
    with open(fname, 'r') as f:
        header = f.readline().rstrip()
        seq = ''
        line = f.readline().rstrip()
        while line != '':
            seq += line
            line = f.readline().rstrip()

    return seq

def bases(seq):
    """
    This function returns the base composition and size of a DNA string.
    PARAMETERS:
        seq - DNA sequence string
    RETURN:
        dict - dictionary of frequencies of each base
        total - DNA sequence total bases dictionary
        perc - dictionary of percentage sof each base in DNA sequence
    """
    a = 0
    c = 0
    t = 0
    g = 0
    count= 0

    for ch in seq:
        if ch == 'A' or ch == 'a':
            a = a + 1
        elif ch == 'C' or ch == 'c':
            c = c + 1
        elif ch == 'T' or ch == 't':
            t = t + 1
        elif ch == 'G' or ch == 'g':
            g = g + 1
        count = count + 1
    pa = (float(a)/len(seq))*100
    pc = (float(c)/len(seq))*100
    pt = (float(t)/len(seq))*100
    pg = (float(g)/len(seq))*100
    dict = {'A': a, 'C': c, 'T': t, 'G': g}
    total = {'Total': count}
    perc = {'Percentage of A': pa, 'Percentage of C': pc, 'Percentage of T': pt, 'Percentage of G': pg}
    return dict, total, perc

def localizedBaseContent(seq, window_size, interval):
    """
    Given a DNA sequence,the function plots a graph of
    percentages of base compositions for a subsequence of size, window_size.
    PARAMATERS:
        seq - DNA sequence string
        window_size - the length DNA sequence being analysed
        interval - the amount by which we move the window_size
    RETURN:
        None
    """

    a_graph = []
    c_graph = []
    t_graph = []
    g_graph = []
    intz = []
    i = 0
    while (i < len(seq) - interval):
        subSeq = seq[i:i+window_size]
        local_base_dict, local_base_total, local_base_percentages = bases(subSeq)

        a_percent = float(local_base_dict['A'])/float(local_base_total['Total'])
        c_percent = float(local_base_dict['C'])/float(local_base_total['Total'])
        t_percent = float(local_base_dict['T'])/float(local_base_total['Total'])
        g_percent = float(local_base_dict['G'])/float(local_base_total['Total'])

        a_graph.append(a_percent)
        c_graph.append(c_percent)
        t_graph.append(t_percent)
        g_graph.append(g_percent)
        intz.append(i)
        i += interval


    plt.figure(figsize=(9,3))
    plt.plot(intz, a_graph, label= 'a')
    plt.plot(intz, c_graph, label= 'c')
    plt.plot(intz, t_graph, label= 't')
    plt.plot(intz, g_graph, label= 'g')
    plt.suptitle('Localized Base Content: ' + 'Window Size = ' + str(window_size) + 'Interval = ' + str(interval))
    plt.legend()
    plt.show()

def skew(seq):
    """
    This function plots the skew diagram for the DNA sequence,
    and returns the position for the minimum skew value.
    PARAMATERS:
        seq - DNA sequence string
    RETURN:
        minimum skew value
    """

    skew_var = 0
    skew_list = []
    for ch in seq:
        if ch == 'C' or ch == 'c':
            skew_var -= 1
        elif ch == 'G' or ch == 'g':
            skew_var += 1
        skew_list.append(skew_var)

    plt.figure(figsize=(9,3))
    plt.plot(range(len(seq)), skew_list)
    plt.suptitle('Skew Diagram')
    plt.show()

    return [str(ch) for ch, j in enumerate(skew_list) if j == min(skew_list)]

def freqkmers(seq, k):
    """
    This function returns the most frequent k-mers in a given DNA sequence
    and a k value.
    PARMETERS:
        seq - DNA sequence string
    RETURN:
        mostFrequentKmers - set of most frequent k-mers in seq
    """
    kmers = {}
    n = len(seq) - k+1
    for i in range(n):
        x = seq[i:i + k]
        if x not in kmers:
            kmers[x] = 0
        kmers[x] +=1
        #kmers.append(x)
    maxF = max(kmers.values())
    mostFrequentKmers = set()
    for x in kmers:
        if kmers[x] == maxF:
            mostFrequentKmers.add(x)
    return mostFrequentKmers

def complement(seq):
    """
    Given a DNA sequence, this function returns the complement
    of the sequence.
    PAREMETERS:
        seq - DNA sequence string
    RETURN:
        comp - complement of DNA sequence
    """
    comp = {'A':'T','T':'A','C':'G','G':'C'}[seq]
    return comp

def revcomplement(seq):
    """
    Given a DNA sequence, this function returns the reverse
    complement of the sequence.
    PAREMETERS:
        seq - DNA sequence string
    RETURN:
        revstr - reverse complement of DNA sequence
    """
    revstr = ""
    for n in seq:
        if n == 'a' or n == 'A':
            comp = 'T'
        elif n == 't' or n == 'T':
            comp = 'A'
        elif n == 'g' or n == 'G':
            comp = 'C'
        elif n == 'c' or n == 'C':
            comp = 'G'
        revstr = revstr + comp
    revstr = revstr[::-1]
    return revstr

def hamming(a, b):
    """
    Given a two values, this function calculates the distance
    between them.
    PAREMETERS:
        a - value 1
        b - value 2
    RETURN:
        c - distance between a & b
    """
    c = 0
    for (x, y) in zip(a, b):
        if x != y:
            c += 1
    return c


def neighbors(dna, d):
    """
    Given a DNA sequence, this function returns the neighbors
    of the sequence with up to d mismatches.
    PAREMETERS:
        seq - DNA sequence string
        d - number of mismatches
    RETURN:
        neighbor - complement of DNA sequence
    """
    if d == 0:
        return dna
    if len(dna) == 1:
        return ["A", "C", "G", "T"]
    neighbor = []
    suffixneighbors = neighbors(dna[1:], d)
    for xy in suffixneighbors:
        if hamming(dna[1:], xy) < d:
            for x in ["A", "C", "G", "T"]:
                neighbor.append(x+xy)
        else:
            neighbor.append(dna[0]+xy)
    return neighbor


def mismatchesandrevcomplements(seq, k, d):
    """
    This function the k-mers that maximize the sum Countd(dna, pattern)
    + Countd(dna, pattern) where pattern is the k-mer and pattern
    is its reverse complement.
    PARMETERS:
        seq - DNA sequence string
        k - length of k-mer/pattern
        d - number of mismatches
    RETURN:
        final - the most frequent k-mers with up to d mismatches
        in a DNA sequence and its reverse complement
    """
    kmers = []
    neighborz = set()
    final = []

    for i in range(len(seq) - k+1):
        kmers.append(seq[i: i+k])
    for i in range(len(seq) - k+1):
        kmers.append(revcomplement(seq[i: i+k]))
    for kmer in kmers:
        neighborz.update(set(neighbors(kmer, d)))
    maxfreq = 0
    for i in neighborz:
        freq = 0
        for km in kmers:
            if hamming(i, km) <= d:
                freq += 1
        if maxfreq < freq:
            maxfreq = freq
            final = [i]
        elif maxfreq == freq:
            final.append(i)
    return final

def main():
    """
    This function calls each function to test it out.
    """
    seq = readFile('sequence.fasta') # question 2
    #base_dict, base_total, base_percentages = bases(seq) # question 3
    #print(base_dict, base_total, base_percentages)
    #localizedBaseContent(seq, 20000, 100) #question 4 a - 20000
    #localizedBaseContent(seq, 90000, 100000) #question 4 b - 90000
    #s = skew(seq) #question 5
    #print(s)
    #print(freqkmers(seq[2562296:2562296+500], 9)) #question 6
    #print(mismatchesandrevcomplements(seq[2562296:2562296+500], 9, 1)) #question 7

main()
