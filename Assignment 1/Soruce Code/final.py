import math
import numpy as np

#==================== GLOBAL VARIABLES
stru_dict = {'H': 0, 'E': 1, 'C': 2}
stru_list = ['H', 'E', 'C']

amino_dict = {'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4, 'GLU': 5, 'GLN': 6, 'GLY': 7, 'HIS': 8,
           'ILE': 9, 'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14, 'SER': 15, 'THR': 16,
           'TRP': 17, 'TYR': 18, 'VAL': 19}

amino_list = np.array(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU','LYS',
              'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'])

amino_total_count = {}

h_cnt = 0
b_cnt = 0
c_cnt = 0


def count(txt_file):
    global h_cnt, b_cnt, c_cnt, amino_total_count
    #new_list = np.empty(shape=(1,5), dtype=str )
    new_list = []
    amino_total_count = {'ALA': 0, 'ARG': 0, 'ASN': 0, 'ASP': 0, 'CYS': 0, 'GLU': 0, 'GLN': 0, 'GLY': 0, 'HIS': 0,
           'ILE': 0, 'LEU': 0, 'LYS': 0, 'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0,
           'TRP': 0, 'TYR': 0, 'VAL': 0}
    h_cnt = 0
    b_cnt = 0
    c_cnt = 0

    txt_file = np.char.replace(txt_file, 'Helix', 'H')
    txt_file = np.char.replace(txt_file, 'Beta', 'E')
    txt_file = np.char.replace(txt_file, 'Other', 'C')
    txt_file = np.char.replace(txt_file, 'Coil', 'C')

    for i in range(len(txt_file)):  
        if txt_file[i][3] in amino_list:
            c = amino_total_count.get(txt_file[i][3])
            c += 1
            amino_total_count.update({txt_file[i][3]: c})

            if txt_file[i][4] == 'H':
                h_cnt += 1
            
            elif txt_file[i][4] == 'E':
                b_cnt += 1
            
            elif txt_file[i][4] == 'C':
                c_cnt += 1

            new_list.append(txt_file[i])
            #new_list = np.append(new_list, [txt_file[i]] , axis=0)
        
    
    return new_list

def separate_protein_class(dataset):
   
    dataset = np.array(dataset, dtype=str)
    _, indices = np.unique(dataset[:, 0], return_index=True)
    indices = np.sort(indices)
    new_list = np.split(dataset, indices)[1:]

    return new_list

# Functions For Self information Table 
def self_create_table():

    table = [[0 for j in range(3)] for i in range(20)]  # 20 proteins are added in the list
    return table

def self_setter(table, residue, structure, count):
    table[residue][structure] = count

def self_getter(table, residue, structure):
    return table[residue][structure]

def self_remove_row(table, stable):     
    for row in range(20):
        for col in range(3):
            table[row][col] -= stable.table[row][col]  

def self_add_row(table, stable):     
    for row in range(20):
        for col in range(3):
            table[row][col] += stable.table[row][col]

#Functions For Residual Information Table
def residual_create_table():    
    table = [[[[0 for n in range(3)] for m in range(20)] for j in range(18)] for i in range(20)]
    return table

def residual_setter(table, j_residue, m, m_residue, j_structure, count):
    table[j_residue][m][m_residue][j_structure] = count

def residual_getter(table, j_residue, m, m_residue, j_structure):
    return table[j_residue][m][m_residue][j_structure]
    
def residual_remove_row(table, ptable):   
    for i1 in range(20):
        for i2 in range(18):
            for i3 in range(20):
                for i4 in range(3):
                    table[i1][i2][i3][i4] -= ptable.table[i1][i2][i3][i4]
                        
def residual_add_row(table, ptable):  
    for i1 in range(20):
        for i2 in range(18):
            for i3 in range(20):
                for i4 in range(3):
                    table[i1][i2][i3][i4] += ptable.table[i1][i2][i3][i4]
    
def create_table(protein_list, self_table, residual_table):
    for j in range(len(protein_list)):
        self_name = protein_list[j][-2]  # the residue name of j
        self_stru = protein_list[j][-1]  # the secondary structure of j
        # update the self information table
        count = self_getter(self_table, amino_dict[self_name], stru_dict[self_stru])
        self_setter(self_table, amino_dict[self_name], stru_dict[self_stru], count + 1)
        # update the pair information table
        sequence0 = protein_list[j][-3]        # the residue sequence code at the location of j
        for m in range(-8, 9, 1):
            if m == 0:  # except for base residue
                continue
            if j + m < 0 or j + m > len(protein_list) - 1:  # avoid index out of range
                continue
            sequence1 = (protein_list[j + m])[-3]       # the residue sequence code at the location of j+m
            m0 = int(sequence1) - int(sequence0)        # the relative location of sequence code
            if m0 not in range(-8, 9, 1):    # avoid sequence breaked
                continue
            residual_name = (protein_list[j + m])[-2]
            # update the pair information table
            count = residual_getter(residual_table, amino_dict[self_name], m0 + 8, amino_dict[residual_name], stru_dict[self_stru])
            residual_setter(residual_table, amino_dict[self_name], m0 + 8, amino_dict[residual_name], stru_dict[self_stru], count + 1)

def generate(aGroup):
    self_table = self_create_table()
    residual_table = residual_create_table()
    
    for protein in aGroup:
        create_table(protein, self_table, residual_table)
    return self_table, residual_table

def compute_self_information(self_table, residue_name):
    I = []      # use to store the I of Helix, Sheet and Coil
    for stru in stru_list:      
        f_sr = self_getter(self_table, amino_dict[residue_name], stru_dict[stru])     
        #         self.table[amino_dict[residue]][stru_dict[structure]]
        f_nsr = 0       
        f_ns = 0        
        f_s = 0         
        for s in stru_list:
            if s != stru:
                f_nsr = f_nsr + self_getter(self_table, amino_dict[residue_name], stru_dict[s])
        for residue in amino_list:
            for s in stru_list:
                if s != stru:
                    f_ns = f_ns + self_getter(self_table, amino_dict[residue], stru_dict[s])
                else:
                    f_s = f_s + self_getter(self_table, amino_dict[residue], stru_dict[s])
        if f_nsr == 0:
            i = math.inf
        else:
            i = math.log(f_sr/f_nsr) + math.log(f_ns/f_s)
        I.append(i)
    return I

def predict(protein, j, self_table, residual_table):    
    self_name = (protein[j])[-2]
    I = compute_self_information(self_table, self_name)  # compute the self-information propensity
    # compute the pair information propensity
    sequence0 = protein[j][2]
    for a in range(3):
        stru = stru_list[a]
        for m in range(-8, 9, 1):
            if m == 0:
                continue
            if j + m < 0 or j + m > len(protein) - 1:  # avoid index out of range
                continue
            sequence1 = (protein[j + m])[2]
            m0 = int(sequence1) - int(sequence0)  # the minus of two residues sequence code
            if m0 not in range(-8, 9):
                continue
            self_name = (protein[j])[-2]
            residual_name = (protein[j + m])[-2]
            number0 = residual_getter(residual_table, amino_dict[self_name], m0 + 8, amino_dict[residual_name], stru_dict[stru])
            number1 = 0
            number2 = 0
            number3 = self_getter(self_table, amino_dict[self_name], stru_dict[stru])
            for s in stru_list:
                if s != stru:
                    number1 = number1 + residual_getter(residual_table, amino_dict[self_name], m0 + 8, amino_dict[residual_name], stru_dict[s])
                    number2 = number2 + self_getter(self_table, amino_dict[residual_name], stru_dict[s])
            if number1 == 0:
                I[a] = math.inf
                break
            try:
                I[a] = I[a] + math.log(number0 / number1) + math.log(number2 / number3)
            except ValueError:
                I[a] = -math.inf
                break
    # check which one is the biggest
    max_value = max(I)
    for i in range(3):
        if I[i] == max_value:
            i = i
            break
    return stru_list[i]

# Predict Both Datasets
def predict_dataset(nested_list, result_list, data_self_table, data_residual_table):
    for i in range(len(nested_list)):
        self_table = self_create_table()
        residual_table = residual_create_table()
        for j in range(len(nested_list)):
            if i == j:
                protein = nested_list[i]
                create_table(protein, self_table, residual_table)
                break
        pred = []
        for j in range(len(protein)):
            pred.append(predict(protein, j, data_self_table, data_residual_table))
        result_list.append(pred)

def score(prediction, facts):
    scores = []
    # traverse both the prediction and actual structure and compare them
    for pChain, fChain in zip(prediction, facts):
        correct = 0
        MCC_table = [[0 for j in range(4)] for i in range(3)]
        for p, f in zip(pChain, fChain):
            if p == f[-1]:
                correct += 1
            for i in range(3):
                residue = stru_list[i]
                if f[-1] == residue:
                    if f[-1] == p:  # TP
                        MCC_table[i][0] += 1
                    else:   # FN
                        MCC_table[i][2] += 1
                else:
                    if p == residue:    # FP
                        MCC_table[i][1] += 1
                    else:
                        MCC_table[i][3] += 1  # TN
        MCC = []
        for table in MCC_table:
            TP = table[0]
            FP = table[1]
            FN = table[2]
            TN = table[3]
            try:
                mcc = (TP*TN - FP*FN)/math.sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
            except ZeroDivisionError:
                mcc = (TP*TN - FP*FN)
            MCC.append(mcc)
        Q3 = correct/len(pChain)
        scores.append([Q3, MCC])
    return scores

# predict the protein family
# it is associated with the MCC.Helix ,MCC.Sheet and MCC.Coil

# prediction for chain
def pFamily(scores):
    def criterion(score):
        # the chain is mainly composed by Alpha
        if abs(score[1]) <= 0.03:
            return 'Alpha'
        # the chain is mainly composed by Beta
        elif abs(score[0]) <= 0.03:
            return 'Beta'
        # the chain is mainly composed by alternating alpha helices and beta seets
        else:
            return 'Alpha/beta'
    
    family = []
    for score in scores:
        family.append(criterion(score[-1]))
    return family

# accuracy of protein family prediction
def accuracy(pFamily, facts):
    count = 0
    for i, j in zip(pFamily, facts):
        if i == j[-1]:
            count += 1
    return count/len(facts)


def main():
    stride_list = np.genfromtxt('/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/stride_info.txt', delimiter="\t", dtype=str)
    dssp_list = np.genfromtxt('/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/dssp_info.txt', delimiter="\t", dtype=str)
    cath_list = np.genfromtxt('/Users/aneruthmohanasundaram/Documents/GitHub/Amino-Acid-Prediction/Assignment 1/Datasets/cath_info.txt', delimiter="\t", dtype=str)

    stride_list = np.char.strip(stride_list)
    dssp_list = np.char.strip(dssp_list)
    cath_list = np.char.strip(cath_list)

    print('\nThe division of structure type in the STRIDE dataset is as follows: ')
    stride_list = count(stride_list)
    print('The Helix count is: ' + str(h_cnt),'\nThe Beta count is:' + str(b_cnt),'\nThe Coil count is:' + str(c_cnt),'\n',amino_total_count)
    
    print('\nThe division of structure type in the DSSP dataset is as follows: ')
    dssp_list = count(dssp_list)
    print('The Helix count is: ' + str(h_cnt),'\nThe Beta count is:' + str(b_cnt),'\nThe Coil count is:' + str(c_cnt),'\n',amino_total_count)

    stride_nest = separate_protein_class(stride_list)
    dssp_nest = separate_protein_class(dssp_list)   

    self_table_stride, residual_table_stride = generate(stride_nest)
    self_table_dssp, residual_table_dssp = generate(dssp_nest)

    pred_stride, pred_dssp = [], []
    stride_result, dssp_result = [], [[]]
    
    predict_dataset(stride_nest, pred_stride, self_table_stride, residual_table_stride)
    predict_dataset(dssp_nest, pred_dssp, self_table_dssp, residual_table_dssp)

    scores_stride = score(pred_stride, stride_nest)
    scores_dssp = score(pred_dssp, dssp_nest)

    count1 = 0
    count2 = 0

    for el1, el2 in zip(scores_stride, scores_dssp):
        if el1[0] >= 0.5:
            count1 += 1
        if el2[0] >= 0.5:
            count2 += 1

    print('\nThe overall Q3 for stride without leaving out protein in the dataset is', 100*count1/len(scores_stride), '%.')
    print('The overall Q3 for dssp without leaving out protein in the dataset is', 100*count2/len(scores_dssp), '%.')

    predf_stride = pFamily(scores_stride)
    predf_dssp = pFamily(scores_dssp)

    print("\nThe accuracy of prediction for protein family in stride dataset is ", accuracy(predf_stride, cath_list))
    print("The accuracy of prediction for protein family in dssp dataset is  ",  accuracy(predf_dssp, cath_list))# predict the protein family   


if __name__ == "__main__":
    main()