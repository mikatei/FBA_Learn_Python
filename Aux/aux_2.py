# Extra functions for FBA

import os


home_path = os.path.dirname(os.path.abspath(__file__))[:-3]
#print("PATH:")
#print(home_path)


def find_index_of_compound(stoichiomatrix_array, compound_name):

    for i in range(len(stoichiomatrix_array)):
        if stoichiomatrix_array[i][0] == compound_name:
            return i

    return -1

#This function takes a file like the ones in the example folder and 
# Creates a d2 list (table list)
def get_rxn_list_d2_example(file_name):
    f = open(file_name, 'r')
    X = f.readlines()
    f.close()

    rxn_list_d2 = []

    for i in range(len(X)):
        rxn = X[i]
        rxn = rxn.replace('\n','')
        y = rxn.split(':')
        #We need to return a list with [[rxn string, rxn name], [rxnstr,rxnname],...]
        if len(y)>1:
            z = [y[1],y[0]]
            rxn_list_d2.append(z)

    

    return rxn_list_d2


#This function gets you upper and lower bounds from the list of reactions.
def give_upper_lower_bounds_list_d2(parsed_rxn_list_d4, bounds_value):
    
    bounds = [[],[]]
    for i in range(len(parsed_rxn_list_d4)):
        if parsed_rxn_list_d4[i][2] == '<=>':
            bounds[0].append(-1*bounds_value)
            bounds[1].append(bounds_value)
        elif parsed_rxn_list_d4[i][2] == '=>':
            bounds[0].append(0)
            bounds[1].append(bounds_value)
        elif parsed_rxn_list_d4[i][2] == '<=':
            bounds[0].append(-1*bounds_value)
            bounds[1].append(0)

        else:
            raise ValueError("Direction of reaction none of '<=>','=>', '<='")
    return bounds


def get_filenames():
    filenames = []

    txtbook_file_name= home_path + 'Examples/textbook_94_example.txt'

    paper_file_name = home_path + 'Examples/paper_example.txt'

    simple_1_file_name = home_path + 'Examples/simple_example.txt'

    simple_2_file_name = home_path + 'Examples/2_simple_example.txt'

    glycolysis_file_name = home_path + 'Examples/Glycolysis_example.txt'


    penthose_phosphate = home_path + 'Examples/penthose-phosphate.txt'

    pp2 =  home_path + 'Examples/penthose-phosphate_exp.txt'

    magic_c = home_path + 'Examples/magic_compound.txt'

    filenames = [paper_file_name, txtbook_file_name, simple_1_file_name, simple_2_file_name, glycolysis_file_name, penthose_phosphate, pp2, magic_c]

    return filenames

