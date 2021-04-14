from math import *
import pandas as pd

# Reading datasets
cath = pd.read_csv("/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/cath_info.txt",delimiter = "\t",header=None, names=["PDB_CODE", "CHAIN_CODE", "Prediction"])

dssp = pd.read_csv("/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/dssp_info.txt", delimiter = "\t",header=None, names=["PDB_CODE", "CHAIN_CODE", "SEQ_POS","AMMINO_ACID","Prediction"])

stride = pd.read_csv("/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/stride_info.txt",delimiter = "\t",header=None, names=["PDB_CODE", "CHAIN_CODE", "SEQ_POS","AMMINO_ACID","Prediction"])

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

amino_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H','ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T','TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

amino_count = {'ALA': 0, 'ARG': 0, 'ASN': 0, 'ASP': 0, 'CYS': 0, 'GLU': 0, 'GLN': 0, 'GLY': 0, 'HIS': 0,
           'ILE': 0, 'LEU': 0, 'LYS': 0, 'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0,
           'TRP': 0, 'TYR': 0, 'VAL': 0}

# Stride dataset
stride.replace('Other','C',inplace=True)
stride.replace('Beta','B',inplace=True)
stride.replace('Helix','H',inplace=True)

# Dssp dataset
dssp.replace('Other','C',inplace=True)
dssp.replace('Coil','C',inplace=True)
dssp.replace('Beta','B',inplace=True)
dssp.replace('Helix','H',inplace=True)

# Change "aList" to desired pandas columns. 
def structure_type(dataset):
    blist = []
    for row in dataset:
        if row in amino_dict.keys():
            c = amino_count.get(row)
            c += 1
            amino_count.update({row: c})
            blist.append(row)
    print(amino_count)
    return blist

def count_values(dataset):
    h_count = dataset.query('Prediction == "H"').Prediction.count()
    b_count = dataset.query('Prediction == "B"').Prediction.count()
    c_count = dataset.query('Prediction == "C"').Prediction.count()
    fs = {} 
    fs['Helix'] = h_count
    fs['Beta'] = b_count
    fs['Coil'] = c_count
    return fs

print('The division of structure type in the STRIDE dataset is as follows: ')
stride_list = structure_type(list(stride['AMMINO_ACID']))
print()
print('The division of structure type in the DSSP dataset is as follows: ')
dssp_list = structure_type(list(dssp['AMMINO_ACID']))
protien_family = list(set(dssp['PDB_CODE']))

print()
# The total count of structure preset in each dataset as follows 
print('Count for STRIDE data set')
print(count_values(stride))
print()
print('Count for DSSP data set')
print(count_values(dssp))

# Change the acid value 
dssp.AMMINO_ACID.replace(amino_dict.keys(),amino_dict.values(),inplace=True)
stride.AMMINO_ACID.replace(amino_dict.keys(),amino_dict.values(),inplace=True)

# We define this function because to iterate over each data row in pandas consumes lots of time so we convert it to a mested list
def fileRead(file_name):
    return [list(i) for i in file_name.values]

# Cerate a self info function
def self_info(dataset,pdb):
    self_info_val = {}
    fs = count_values(dataset)
    helix_fs = fs['Helix']
    sheet_fs = fs['Beta']
    coil_fs = fs['Coil']

    for idx,row in dataset.iterrows():
        if row[0] == pdb:
            curr = row[3]
            if row[3] in list(amino_dict.values()):
                helix_fsr,sheet_fsr,coil_fsr = 0,0,0
            
            for idx,other in dataset.iterrows():
                if other[0] != pdb:
                    if other[3] == curr:
                        if other[4] == 'H':
                            helix_fsr += 1
                        if other[4] == 'B':
                            sheet_fsr += 1
                        if other[4] == 'C':
                            coil_fsr += 1
            helix_i = log(helix_fsr / ((helix_fsr + sheet_fsr + coil_fsr) - helix_fsr)) + log((sheet_fs + coil_fs) / helix_fs)
            sheet_i = log(sheet_fsr / ((helix_fsr + sheet_fsr + coil_fsr) - sheet_fsr)) + log((helix_fs + coil_fs) / sheet_fs)
            coil_i = log(coil_fsr / ((helix_fsr + sheet_fsr + coil_fsr) - coil_fsr)) + log((helix_fs + sheet_fs) / coil_fs)
            
            self_info_val[curr] = (helix_i, sheet_i, coil_i)
    return self_info_val

def _calc_pair_info(dataset,pdb):
    '''
    Made to calculate acid pair information.
    Returns dictionary {acid:{position:{acid:(helix_info, sheet_info, coil_info)}}}
    '''
    # acids_list = self.acids[:]
    pair_info = {}
    # rows,cols = dataset.shape
    
    fs = count_values(stride)
    helix_fs = fs['Helix']
    sheet_fs = fs['Beta']
    coil_fs = fs['Coil']
    
    for current_acid in list(amino_dict.values()):
        '''loop over 20 natural amino acids'''
        acid_dict = {}
        
        for m in [-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8]:
            '''loop over 16 neighbors'''
            position_dict = {}
            
            
            for pair_acid in list(amino_dict.values()):
                '''loop over 20 natural amino acids'''
                    
                helix_fsr = 0
                sheet_fsr = 0
                coil_fsr = 0
                
                for i in range(len(dataset) - 8):
                    '''loop over all rest recidues'''
                    if pair_acid == dataset[i+m][3] and dataset[i+m][0] == dataset[i][0] and dataset[i][0] != pdb and dataset[i][3] == current_acid:
                        if dataset[i][4] == 'H':
                            helix_fsr += 1
                        if dataset[i][4] == 'B':
                            sheet_fsr += 1
                        if dataset[i][4] == 'C':
                            coil_fsr += 1

                if ((helix_fsr + sheet_fsr + coil_fsr) - helix_fsr) != 0:
                    helix_log_arg = helix_fsr / ((helix_fsr + sheet_fsr + coil_fsr) - helix_fsr)
                else:
                    helix_log_arg = 1

                if ((helix_fsr + sheet_fsr + coil_fsr) - sheet_fsr) != 0:
                    sheet_log_arg = sheet_fsr / ((helix_fsr + sheet_fsr + coil_fsr) - sheet_fsr)
                else:
                    sheet_log_arg = 1

                if ((helix_fsr + sheet_fsr + coil_fsr) - coil_fsr) != 0:
                    coil_log_arg = coil_fsr / ((helix_fsr + sheet_fsr + coil_fsr) - coil_fsr)
                else:
                    coil_log_arg = 1

                helix_i = log(helix_log_arg if helix_log_arg > 0 else 1) + log((sheet_fs + coil_fs) / helix_fs)
                sheet_i = log(sheet_log_arg if sheet_log_arg > 0 else 1) + log((helix_fs + coil_fs) / sheet_fs)
                coil_i = log(coil_log_arg if coil_log_arg > 0 else 1) + log((helix_fs + sheet_fs) / coil_fs)

                position_dict[pair_acid] = (helix_i, sheet_i, coil_i)
            
            acid_dict[m] = position_dict
        pair_info[current_acid] = acid_dict
    
    return pair_info

def sec_struc_prediction(dataset,PDB_code):
    '''
    Returns a secondary stucture prediction of the protein as a list:
    [PDB_code, the secondary structure prediction, Q3]
    '''
    selfinfo = self_info(dataset,PDB_code)
    pairinfo = _calc_pair_info(dataset,PDB_code)
    protein_list = [entry for entry in self.dataset if entry[0] == PDB_code]
    protein = ''
    structure_real = ''
    prediction = [PDB_code]
    structure_prediction = ''
    
    '''get the real protein secondary structure as a string'''
    for acid in protein_list:
        protein += acid[3]
        
        if acid[4] == 'H':
            structure_real += 'H'
        elif acid[4] == 'B':
            structure_real += 'E'
        elif acid[4] == 'C':
            structure_real += 'C'
    
    # print(protein, structure_real)
    
    '''get the protein secondary stucture prediction as a string'''
    for i in range(len(protein)):
        helix = [selfinfo[protein[i]][0], 'H']
        sheet = [selfinfo[protein[i]][1], 'E']
        coil = [selfinfo[protein[i]][2], 'C']
        
        for m in [-8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8]:
            if i+m > 0 and i+m < len(protein):
                helix[0] += pairinfo[protein[i]][m][protein[i+m]][0]
                sheet[0] += pairinfo[protein[i]][m][protein[i+m]][1]
                coil[0] += pairinfo[protein[i]][m][protein[i+m]][2]
        
        if max(helix[0], sheet[0], coil[0]) == helix[0]:
            structure_prediction += helix[1]
        if max(helix[0], sheet[0], coil[0]) == sheet[0]:
            structure_prediction += sheet[1]
        if max(helix[0], sheet[0], coil[0]) == coil[0]:
            structure_prediction += coil[1]
    
    prediction.append(structure_prediction)
    
    '''Q3 calculation'''
    recidues_predicted = 0
    for i in range(len(structure_real)):
        if structure_real[i] == structure_prediction[i]:
            recidues_predicted += 1
    
    q3 = recidues_predicted / len(structure_real)
    prediction.append(q3)
    
    return prediction

sec_struc_prediction(stride,protien_family)
