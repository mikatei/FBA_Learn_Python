#Adding modularity with small functions.

from Aux.parsing import tsv_to_d2_list
from Aux.aux_2 import find_index_of_compound
import logging



#Finding direction of reaction
#Returns direction of reaction: '<=', '<=>', or '=>'
def reaction_direction(reaction_str):

    #t stands for direction of reaction, can either be '<=', '=>', or '<=>'
    if reaction_str.count('<') + reaction_str.count('>') == 2:
        t = '<=>'
    elif reaction_str.count('<') == 1:
        t = '<='
    elif reaction_str.count('>') == 1:
        t = '=>'
    else:
        t = '0'
        logging.debug("ERROR: incorrectly formatted reaction. Must include '<' or '>' or both.")
        logging.debug(reaction_str)

    return t


def turn_reaction_half_into_list(half_rxn_str):
        output_list = []
        compounds = half_rxn_str.split(' + ')
        for i in range(len(compounds)):
            compound = compounds[i]
            #find index of '(' and ')' then add int(in between)
            loc_1 = compound.find('(')
            loc_2 = compound.find(')')
            #TD Catch errors in loc_1 and loc_2
            c = compound[loc_1 + 1:loc_2]
            if c != '':
                try:
                    num_compound = float(c)
                except ValueError:
                    print("PROBLEM VALUE: " + c + ':' + str(len(c)))
                    print(half_rxn_str)
                    num_compound = "?"

                if num_compound != "?":
                    compound_name = compound[loc_2 + 1:].replace(' ', '')
                    output_list.append([num_compound, compound_name])
            else:
                logging.critical('compound not found')
        return output_list


def extract_compounds_from_rxn_list_d3(rxn_list_d3):
    compounds = []
    substrates = rxn_list_d3[0]
    products = rxn_list_d3[1]

    #substrates looks like [[compound_a],[compound_b], ... ]
    #where each compound looks like [num, compound name], so compound name is in index 1
    for compound in substrates:
        if len(compound) > 0:
            compounds.append(compound[1])
    for compound in products:
        if len(compound) > 0:
            compounds.append(compound[1])

    return compounds


def extract_all_compounds_from_parsed_rxn_list(parsed_rxn_list_d4):

    #At first, compounds will be a set, in order to avoid duplicating elements
    compounds = set()
    for parsed_rxn in parsed_rxn_list_d4:
        new_compounds = extract_compounds_from_rxn_list_d3(parsed_rxn)
        compounds.update(new_compounds)

    final_compounds = list(compounds)
    #logging.debug(len(final_compounds))
    final_compounds.sort()
    if "" in final_compounds:
        final_compounds.remove("")
    return final_compounds
    

def get_rxn_list_d2_from_file(file_name):

    total_table = tsv_to_d2_list(file_name)
    
    #0th column is name; 9th column is the equations, starting at row 1 (row 0 is 'equations')

    rxn_list_d2 = []

    #TEST:

    for i in range(1,10):
        crnt_row = total_table[i]
        rxn_list_d2.append([crnt_row[8],crnt_row[0]])
    '''
    for i in range(1,len(total_table)-1):
        crnt_row = total_table[i]
        rxn_list_d2.append([crnt_row[8],crnt_row[0]])
    '''
    return rxn_list_d2


# Input to this is the parsed reaction list [[substrates],[products],type,name] and the 
# Stoichiomatrix dictionary, where per each compound there is a list [0,0,0,0,0,...] with a 0
# per each reaction.
# Since index i runs through the reactions, the ith row represents it as a reaction.
# UNKNOWN: How to add compound numbers when it comes to reactions <= or <=>
def fill_in_stoichiomatrix_dict(parsed_rxn_list_d4, stoichiomatrix):

    success_var = 0
    for i in range(len(parsed_rxn_list_d4)):
        crnt_rxn = parsed_rxn_list_d4[i]
        crnt_direction = crnt_rxn[2]
        if (crnt_direction == '=>') or (crnt_direction == '<=>'):
            #In this case, the substrates are turning into negatives, and the products positive
            crnt_substrates = crnt_rxn[0]
            for substrate in crnt_substrates:
                if isfloat(substrate[0]):
                    compound_number = float(substrate[0])
                    compound_name = substrate[1]
                    
                    logging.debug(compound_name)
                    if len(compound_name) > 0:
                        new_compound_list = stoichiomatrix[compound_name][:]
                        new_compound_list[i] = new_compound_list[i] - compound_number
                        #logging.debug('A')
                        #logging.debug(new_compound_list)
                        stoichiomatrix.update({compound_name: new_compound_list})
                        #logging.debug('B')
                        #logging.debug(stoichiomatrix[compound_name])
            crnt_products = crnt_rxn[1]
            for product in crnt_products:
                if isfloat(product[0]):
                    compound_number = float(product[0])
                    compound_name = product[1]
                    if len(compound_name) > 0 :
                        new_compound_list = stoichiomatrix[compound_name][:]
                        new_compound_list[i] = new_compound_list[i] + compound_number
                        stoichiomatrix[compound_name] = new_compound_list
        elif crnt_direction == '<=':
            #In this case, the substrates are turning into positives, and the products negative
            crnt_substrates = crnt_rxn[0]
            for substrate in crnt_substrates:
                if isfloat(substrate[0]):
                    compound_number = float(substrate[0])
                    compound_name = substrate[1]
                    if len(compound_name) > 0:
                        new_compound_list = stoichiomatrix[compound_name][:]
                        new_compound_list[i] = new_compound_list[i] + compound_number
                        stoichiomatrix[compound_name] = new_compound_list
            crnt_products = crnt_rxn[1]
            for product in crnt_products:
                if isfloat(product[0]):
                    compound_number = float(product[0])
                    compound_name = product[1]
                    if len(compound_name) > 0:
                        new_compound_list = stoichiomatrix[compound_name][:]
                        new_compound_list[i] = new_compound_list[i] - compound_number
                        stoichiomatrix[compound_name] = new_compound_list
        else:
            logging.debug("ERROR: direction value: '<=>' is not recognized.")
            success_var = 1

    return [stoichiomatrix, success_var]


def fill_in_stoichiomatrix_array(parsed_rxn_list_d4, stoichiomatrix_array):

    success_var = 0
    logging.debug("LEN: " + str(len(parsed_rxn_list_d4)))
    for i in range(1, len(parsed_rxn_list_d4)):
        crnt_rxn = parsed_rxn_list_d4[i]
        crnt_direction = crnt_rxn[2]
        logging.debug(crnt_direction)
        if (crnt_direction == '=>') or (crnt_direction == '<=>'):
            #In this case, the substrates are turning into negatives, and the products positive
            crnt_substrates = crnt_rxn[0]
            for substrate in crnt_substrates:
                if isfloat(substrate[0]):
                    compound_number = float(substrate[0])
                    compound_name = substrate[1]
                    #logging.debug(compound_name)
                    compound_index = find_index_of_compound(stoichiomatrix_array, compound_name)
                    if compound_index != -1:
                        new_compound_list = stoichiomatrix_array[compound_index][:]
                        new_compound_list[i] = float(new_compound_list[i]) - compound_number
                        logging.debug('NCL START')
                        logging.debug(new_compound_list)
                        logging.debug('NCL STOP')
                        stoichiomatrix_array[compound_index] = new_compound_list
            crnt_products = crnt_rxn[1]
            for product in crnt_products:
                if isfloat(product[0]):
                    compound_number = float(product[0])
                    compound_name = product[1]
                    compound_index = find_index_of_compound(stoichiomatrix_array, compound_name)
                    if compound_index != -1:
                        new_compound_list = stoichiomatrix_array[compound_index][:]
                        new_compound_list[i] = float(new_compound_list[i]) + compound_number
                        stoichiomatrix_array[compound_index] = new_compound_list
        elif crnt_direction == '<=':
            #In this case, the substrates are turning into positives, and the products negative
            crnt_substrates = crnt_rxn[0]
            for substrate in crnt_substrates:
                if isfloat(substrate[0]):
                    compound_number = float(substrate[0])
                    compound_name = substrate[1]
                    compound_index = find_index_of_compound(stoichiomatrix_array, compound_name)
                    if compound_index != -1:
                        new_compound_list = stoichiomatrix_array[compound_index][:]
                        new_compound_list[i] = float(new_compound_list[i]) + compound_number
                        stoichiomatrix_array[compound_index] = new_compound_list
            crnt_products = crnt_rxn[1]
            for product in crnt_products:
                if isfloat(product[0]):
                    compound_number = float(product[0])
                    compound_name = product[1]
                    compound_index = find_index_of_compound(stoichiomatrix_array, compound_name)
                    if compound_index != -1:
                        new_compound_list = stoichiomatrix_array[compound_index][:]
                        new_compound_list[i] = float(new_compound_list[i]) + compound_number
                        stoichiomatrix_array[compound_index] = new_compound_list
        else:
            success_var = 1

        return [stoichiomatrix_array, success_var]



def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
