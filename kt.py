import math,pandas as pd

# Reading datasets
cath = pd.read_csv("/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/cath_info.txt",delimiter = "\t",header=None, names=["PDB_CODE", "CHAIN_CODE", "Prediction"])

dssp = pd.read_csv("/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/dssp_info.txt", delimiter = "\t",header=None, names=["PDB_CODE", "CHAIN_CODE", "SEQ_POS","AMMINO_ACID","Prediction"])

stride = pd.read_csv("/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/stride_info.txt",delimiter = "\t",header=None, names=["PDB_CODE", "CHAIN_CODE", "SEQ_POS","AMMINO_ACID","Prediction"])

# Assigning the names 
dssp.name = "DSSP"
stride.name = "STRIDE"

# To count the acids present 
stru_dict = {'H': 0, 'E': 1, 'C': 2}

stru_list = ['H', 'E', 'C']

amino_dict = {'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4, 'GLU': 5, 'GLN': 6, 'GLY': 7, 'HIS': 8,
           'ILE': 9, 'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14, 'SER': 15, 'THR': 16,
           'TRP': 17, 'TYR': 18, 'VAL': 19}

amino_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU','LYS',
              'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

amino_count = {'ALA': 0, 'ARG': 0, 'ASN': 0, 'ASP': 0, 'CYS': 0, 'GLU': 0, 'GLN': 0, 'GLY': 0, 'HIS': 0,
           'ILE': 0, 'LEU': 0, 'LYS': 0, 'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0,
           'TRP': 0, 'TYR': 0, 'VAL': 0}

# Check dataset head
print("Stide Dataset")
print(stride.head())
print()
print("DSSP Dataset")
print(dssp.head())
print()
print("Cath Dataset")
print(cath.head())

# Change "aList" to desired pandas columns. 
def structure_type(dataset):
    blist = []
    for row in dataset:
        if row in amino_list:
            c = amino_count.get(row)
            c += 1
            amino_count.update({row: c})
            # if row == 'Helix':
            #     # row[-1] = 'H'
            #     h_count += 1
            # elif row == 'Beta':
            #     # row[-1] = 'E'
            #     b_count += 1
            # elif row == 'Coil' or row == 'Other':
            #     # row[-1] = 'C'
            #     c_count += 1
            blist.append(row)
    print(amino_count)
    return blist

def count_values(dataset):
    h_count = dataset.query('Prediction == "Helix"').Prediction.count()
    b_count = dataset.query('Prediction == "Beta"').Prediction.count()
    c_count = dataset.query('Prediction == "Other"' or 'Prediction == "Coil"').Prediction.count()
    print(f"{dataset.name} Structures Count:")
    print('The Helix count is: ' + str(h_count))
    print('The Beta count is: ' + str(b_count))
    print('The Coil count is: ' + str(c_count))

# Classifying the amino acid present 
print('The division of structure type in the STRIDE dataset is as follows: ')
stride_list = structure_type(list(stride['AMMINO_ACID']))
amino_count = {'ALA': 0, 'ARG': 0, 'ASN': 0, 'ASP': 0, 'CYS': 0, 'GLU': 0, 'GLN': 0, 'GLY': 0, 'HIS': 0,'ILE': 0, 'LEU': 0, 'LYS': 0,'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0,'TRP': 0, 'TYR': 0, 'VAL': 0}

print()
print('The division of structure type in the DSSP dataset is as follows: ')
dssp_list = structure_type(list(dssp['AMMINO_ACID']))

print()
count_values(stride)
print()
count_values(dssp)