# Defining the path for all our dataset 
cath_path = '/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/cath_info.txt'
dssp_path = '/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/dssp_info.txt'
stride_path = '/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/stride_info.txt'

amino = {'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4, 'GLU': 5, 'GLN': 6, 'GLY': 7, 'HIS': 8,'ILE': 9, 'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14, 'SER': 15, 'THR': 16,'TRP': 17, 'TYR': 18, 'VAL': 19,'UNK': 20}

# To keep a count on amino acid present in our dataset
aminoCount = {'ALA': 0, 'ARG': 0, 'ASN': 0, 'ASP': 0, 'CYS': 0, 'GLU': 0, 'GLN': 0, 'GLY': 0, 'HIS': 0,'ILE': 0, 'LEU': 0, 'LYS': 0, 'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0,'TRP': 0, 'TYR': 0, 'VAL': 0,'UNK': 0}


# Creating our dictionary to have check of type of structure in our dataset such as helix, coil and beta.
stru_dict = {'H': 0, 'E': 1, 'C': 2}
stru_list = ['H', 'E', 'C']

# Since we have a txt file so we need to convert it to a list and choosing list is because they can be iterateable. 
def fileRead(path):
    fileOpen = open(path)
    return [i.split() for i in fileOpen.readlines()]

# This function prints total .... present in our dataset
def structure_type(alist):
    blist = []
    helixCount,betaCount,coilCount = 0,0,0 # Initialisng our helix, beta, coli count as 0 and thus incremented whenever we lookup for it in our dataset. 
    for i in alist: # Traverse through our list where we pass our fileRead function which return "list".
        if i[3] in amino.keys(): # Checks if our thrid index of our list contains Amino acid in our dataset which it compares with our dictionary present.
            c = aminoCount.get(i[3])
            c += 1
            aminoCount.update({i[3]: c}) # If Amino acid matches with our dictionary then we update the values of Amino acid.
            if i[4] == 'Helix':
                i[-1] = 'H'
                helixCount += 1
            elif i[4] == 'Beta':
                i[-1] = 'E'
                betaCount += 1
            elif i[4] == 'Coil' or i[4] == 'Other':
                i[-1] = 'C'
                coilCount += 1
            blist.append(i)
    print(f'The Helix count is: {str(helixCount)}')
    print(f'The Beta count is: {str(betaCount)}')
    print(f'The Coil count is: {str(coilCount)}')
    print(f'Amino Acid count present in our dictionary is: {aminoCount}')
    return blist

# Structure type of Stride dataset path
print('Stride Structure:')
stride_struct = structure_type(fileRead(stride_path))
print('\n')
print('Dssp Structure:')
dssp_struct = structure_type(fileRead(dssp_path))