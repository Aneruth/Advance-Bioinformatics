# Importing packages to visualise
import matplotlib.pyplot as plt

genome_info = '/Users/aneruthmohanasundaram/Documents/GitHub/Advance-Bioinformatics/Assignment 2/Datasets/dm6.fa'
sequence_info = '/Users/aneruthmohanasundaram/Documents/GitHub/Advance-Bioinformatics/Assignment 2/Datasets/10k_reads.fastq'

def readGenome(aFile):
    ''' Funtion to read the fastq file '''
    f = open(aFile)
    aval = ''.join([i.rstrip() for i in f.readlines() if not i[0] == '>'])
    return aval.upper()

def countFrequency(listToPass):
    ''' This functions counts the frequency of each base present in our genome '''
    base_count = {'A':0,'G':0,'C':0,'T':0}
    for i in listToPass:
        if i in base_count.keys():
            base_count[i] += 1
    return base_count

def readSequence(seq):
    f = open(seq,'r')
    sequence,quality = [],[]
    while True:
        f.readline()
        seq = f.readline().rstrip()
        f.readline()
        qual = f.readline().rstrip()
        if len(seq) == 0: break # Exit the while loop
        sequence.append(seq)
        quality.append(qual)
    return sequence,quality

genome = readGenome(genome_info)
asd = genome[:500]
print(asd)
print(countFrequency(asd))