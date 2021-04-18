genome_info = '/Users/aneruthmohanasundaram/Documents/GitHub/Advance-Bioinformatics/Assignment 2/Datasets/dm6.fa'
sequence_info = '/Users/aneruthmohanasundaram/Documents/GitHub/Advance-Bioinformatics/Assignment 2/Datasets/10k_reads.fastq'

def fileRead(aFile):
    ''' Funtion to read the fastq file '''
    f = open(aFile)
    aval = ''.join([i.rstrip() for i in f.readlines() if not i[0] == '>'])
    return aval

genome = fileRead(genome_info)
asd = genome[:500]

def countFrequency(listToPass):
    ''' This functions counts the frequency of each base present in out genome '''
    base_count = {'A':0,'G':0,'C':0,'T':0}
    for i in listToPass:
        if i in base_count.keys():
            base_count[i] += 1
    return base_count
print(asd)
print(countFrequency(asd))