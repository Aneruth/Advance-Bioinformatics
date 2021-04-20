# Importing packages to visualise
import matplotlib.pyplot as plt

genome_info = '/Users/aneruthmohanasundaram/Documents/GitHub/Advance-Bioinformatics/Assignment 2/Datasets/dm6.fa'
sequence_info = '/Users/aneruthmohanasundaram/Documents/GitHub/Advance-Bioinformatics/Assignment 2/Datasets/10k_reads.fastq'

def readGenome(aFile):
    ''' Funtion to read the fastq file. In our fastq file we have a special character and we need to consider our codes after that. '''
    f = open(aFile)
    aval = ''.join([i.rstrip() for i in f.readlines() if not i[0] == '>'])
    return aval.upper() # Since our characters consist of mixed case letters.

def countFrequency(listToPass):
    ''' This functions counts the frequency of each base present in our genome '''
    base_count = {'A':0,'G':0,'C':0,'T':0} # To keep a count number of bases present in our dataset
    for i in listToPass:
        if i in base_count.keys():
            base_count[i] += 1
    return base_count

def readSequence(seq):
    ''' A function to extract all the sequence and quality score from our dataset '''
    f = open(seq,'r') # Reading the file
    sequence,quality = [],[]
    while True:
        f.readline()
        seq = f.readline().rstrip() # Assigning our sequence to a variable
        f.readline()
        qual = f.readline().rstrip() # Fetching our quality 
        if len(seq) == 0: break # Exit the while loop
        sequence.append(seq) # Appending our sequence to our newly created empty list
        quality.append(qual) # Appending our quality score to newly created empty list
    return sequence,quality # Returns the values in tuple format

def qulaityScore(astring):
    return ord(astring) - 33 # Ord function gives the ascii values

def fetchQualityScore(alist):
    blist = [qulaityScore(j) for i in list(alist) for j in i]
    return blist

def suffixArray(alist):
    pass

genome = readGenome(genome_info)
asd = genome[:500] # To read the first 500 genome 
sequence_list = readSequence(sequence_info)[0]
quality_list = readSequence(sequence_info)[1]
qulaity_score = fetchQualityScore(quality_list)

# print(sequence_list[:5])
print(qulaity_score[:10])
# print(countFrequency(asd))

# Visulaising our dataset
keys = countFrequency(asd).keys()
values = countFrequency(asd).values()
plt.figure(figsize=(15,10))
plt.bar(keys,values)
plt.title('Number of Bases present')
plt.show()