# This file contains the functions required to build a stoichiometric matrix.
# Info:
#   To construct the stoichiometric matrix, all compounds used in all reactions are stored
#   in  the columns of the matrix, and the individual cells contain the stoichiometry of each 
#   compound in each reaction: negative values indicate compound is consumed in the reaction.
#   Zero indicating it is not involved in the reaction. Positive indicate produced by reaction.
#   (Most values are zero). Dimensions are number of compounds and number of reactions.
#   Each pathway could have its own stoichiometric matrix.



# Build stoichiometric matrix function takes a list of reactions and compounds and builds a stoichiometric matrix.
# Reactions are notated by: (#) cpdA + (#) cpdB [<]=[>] (#) cpdC + (#) cdpD

from Aux.aux_1 import reaction_direction, turn_reaction_half_into_list, extract_all_compounds_from_parsed_rxn_list, get_rxn_list_d2_from_file, fill_in_stoichiomatrix_dict, fill_in_stoichiomatrix_array
from Aux.aux_2 import get_rxn_list_d2_example
import logging
import numpy as np
from scipy import linalg





#Supposing a reaction is given in the format: 
#    (1)  cpd00012[c0] + (1)  cpd02557[c0] <= (1)  cpd00113[c0] + (1)  cpd02590[c0]
# We can split it by first counting the number of '<' and '>' and then splitting by the appropriate
# intermediate term. 
# Output looks like: [[inputs], [outputs]] where [inputs] = [[compound_a],[compound_b],...]
# And each compound looks like [compound_x] = [#, name]. This means output is d3?
# We also want the name of the reaction


# The input to this will be a list of [reaction_string, reaction_name]
# The output is either a single parsed reaxtion in a list
def parse_reaction_into_compound_numbers(reaction_list_d1, therm_bool):

    parsed_rxns_subset_d4 = []
    reaction_str = reaction_list_d1[0]
    reaction_name = reaction_list_d1[1]
    #t stands for direction of reaction, can either be '<=', '=>', or '<=>'
    t = reaction_direction(reaction_str)
    
    #Splitting the reaction string into substrates and products strings.
    split_1 = reaction_str.split(t)
    substrates_str = split_1[0]
    products_str = split_1[1]

    #Now we pull out each compound and its number
    substrates = turn_reaction_half_into_list(substrates_str)
    products = turn_reaction_half_into_list(products_str)
    
    #Depending on whether we are using thermodynamics or not, you may need to split the bidirectional reactions into
    #    separate components (i.e. <=> turns into both => and <=)
    if therm_bool == True:

        if t == '<=>':
            parsed_rxns_subset_d4.append([substrates,products,'=>',reaction_name + 'a'])
            parsed_rxns_subset_d4.append([products,substrates,'=>',reaction_name + 'b'])
            return parsed_rxns_subset_d4
        else:
            parsed_rxns_subset_d4.append([substrates, products, t, reaction_name])
            return parsed_rxns_subset_d4
    else:
            parsed_rxns_subset_d4.append([substrates, products, t, reaction_name])
            return parsed_rxns_subset_d4




        
#NOTE: in the downloaded model, the last reaction is different, eg 'bio1'- which should not be used.
def list_of_reaction_strings_to_parsed_reaction_list(rxn_string_list_d2, therm_bool):
    parsed_rxn_list_d4 = []
    for i in range(len(rxn_string_list_d2)):
        parsed_rxns_sub_d4 = parse_reaction_into_compound_numbers(rxn_string_list_d2[i], therm_bool)
        for parsed_rxn_list_d3 in parsed_rxns_sub_d4:
            parsed_rxn_list_d4.append(parsed_rxn_list_d3)
    return parsed_rxn_list_d4

def get_rxn_name_list(parsed_rxn_list_d4):
    rxn_name_list = [' ']
    for i in range(len(parsed_rxn_list_d4)):
        rxn_name = parsed_rxn_list_d4[i][3]
        rxn_name_list.append(rxn_name)
    return rxn_name_list



#From the parsed reactions of a genome, we can create the stoichiometric matrix
#The parsed rxn list looks like: [[rxn_1],[rxn_2],[rxn_3],...]. Each rxn looks like:
# [[substrates],[products],direction (<=, =>, <=>), name] 
def create_stoichiometric_matrix(parsed_rxn_list_d4):


    # 1 column per reaction
    num_of_columns = len(parsed_rxn_list_d4)

    #compounds is just a list of strings with compound names
    compounds = extract_all_compounds_from_parsed_rxn_list(parsed_rxn_list_d4)
    num_of_rows = len(compounds)

    #We'll make the stoichiomatrix a dict.
    stoichiomatrix_dict = {}


    #Trying with an array instead
    stoichiomatrix_array = []


    for i in range(num_of_rows):
        compound = compounds[i]
        stoichiomatrix_array.append([compound] + [0]*num_of_columns)
        stoichiomatrix_dict[compound] = [0]*num_of_columns
    
    #Now the stoichiomatrix is initialized. We need to fill in the locations for the reactions.
    stoichiomatrix_dict_updated = fill_in_stoichiomatrix_dict(parsed_rxn_list_d4, stoichiomatrix_dict)[0]

    #We skip the array option for now.
    #stoichiomatrix_array_updated = fill_in_stoichiomatrix_array(parsed_rxn_list_d4, stoichiomatrix_array)[0]

    
    #TEST:
    #print('COMPOUNDS-------------------')
    #Calculating non-zero elements
    c = 0
    comp_count = 0
    max_count = 0
    max_compound = ''
    for i in range(len(compounds)):
        comp_count = 0
        compound = compounds[i]
        X = stoichiomatrix_dict_updated[compound]
        #print(X)
        for j in range(0, len(X)):
            if X[j] != 0:
                c += 1
                comp_count += 1
        if comp_count > max_count:
            max_count = comp_count
            max_compound = compound
       
    
    #TEST 2
    #print(stoichiomatrix_dict_updated['cpd02074[c0]'])
    
    logging.info("Total Elements in Matrix: " + str(len(compounds)*len(parsed_rxn_list_d4)))
    logging.info("Total Reaction Number: " + str(len(parsed_rxn_list_d4)))
    logging.info("Total Compounds: " + str(len(compounds)))
    logging.info("Total nonzero elements: " + str(c))
    logging.info("Percentage: " + str(float(c/(len(compounds)*len(parsed_rxn_list_d4)))))

    logging.info("max count: " + str(max_count))
    logging.info("max compound: " + max_compound)

    #cmpnd_mtrx holds the compound names as well as the stoichiometric values
    cmpnd_mtrx = []
    #whereas mtrx only holds the stoichiometric values
    mtrx = []

    #first, we add the reaction names to the matrix.
    cmpnd_mtrx.append(get_rxn_name_list(parsed_rxn_list_d4))

    for k in stoichiomatrix_dict_updated.keys():
        cmpnd_mtrx.append([k] + stoichiomatrix_dict_updated[k])
        mtrx.append(stoichiomatrix_dict_updated[k])

    S_with_compounds = np.array(cmpnd_mtrx)
    S = np.array(mtrx)


    #return stoichiomatrix_dict_updated
    return [S_with_compounds,S]


#Testing on simple example inputs
def test_sm():
    paper_file_name='/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA/Examples/paper_example.txt'
    file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA/Examples/textbook_94_example.txt'
    rxn_list_d2 = get_rxn_list_d2_example(paper_file_name)
     
    parsed_rxn_list_d4 = list_of_reaction_strings_to_parsed_reaction_list(rxn_list_d2)

    mtrices = create_stoichiometric_matrix(parsed_rxn_list_d4)
    S_w_cmpnds = mtrices[0]
    S = mtrices[1]
    
    print(S_w_cmpnds)
    print(S)

def main():

    test_sm()

#main()
