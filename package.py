# Predict Both Datasets
def predict_dataset(nested_list, result_list, data_self_table, data_residual_table):
    for i in range(len(nested_list)):
        self_table = self_create_table()
        residual_table = residual_create_table()
        for j in range(len(nested_list)):
            if i == j:
                protein_list = nested_list[i]
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
                break
        pred = []
        for j in range(len(protein_list)):
            pred.append(predict(protein_list, j, data_self_table, data_residual_table))
        result_list.append(pred)

def score(prediction, facts):
    scores = []
    # traverse both the prediction and actual structure and compare them
    for pChain, fChain in zip(prediction, facts):
        correct = 0
        MCC_table = np.zeros(shape=(3,4),dtype=int).tolist()
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