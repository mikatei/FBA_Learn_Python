'''
This file holds functions that check the FBA process for imbalances (e.g. A compound is used up in
reactions, but has no way to enter the cell).
'''



#This function uses the Stoichiometric Matrix as an input
# It returns a list of imbalanced_compounds
def get_indices_of_imbalanced_compounds(Mt):
    imbalanced_compounds = []

    for i in range(len(Mt)):
        cmpnd_row = Mt[i]

        #We initially assume the compound is neither produced or used. They should both be true.
        produced = False
        used = False
        for j in range(len(cmpnd_row)):
            if cmpnd_row[j] > 0:
                produced = True
            elif cmpnd_row[j] < 0:
                used = True
        if produced == False or used == False:
            imbalanced_compounds.append(i)

    return imbalanced_compounds




