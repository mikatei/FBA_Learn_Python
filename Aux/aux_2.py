# Extra functions for FBA



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

    paper_file_name='/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA_Learn_Python/Examples/paper_example.txt'
    txtbook_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA_Learn_Python/Examples/textbook_94_example.txt'

    simple_1_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA_Learn_Python/Examples/simple_example.txt'

    simple_2_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA_Learn_Python/Examples/2_simple_example.txt'

    glycolysis_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA_Learn_Python/Examples/Glycolysis_example.txt'


    penthose_phosphate = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA_Learn_Python/Examples/penthose-phosphate.txt'

    pp2 =  '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA_Learn_Python/Examples/penthose-phosphate_exp.txt'

    filenames = [paper_file_name, txtbook_file_name, simple_1_file_name, simple_2_file_name, glycolysis_file_name, penthose_phosphate, pp2]

    return filenames

