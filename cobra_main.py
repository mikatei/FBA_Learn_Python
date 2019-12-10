#python




import cobra.test
import cobra
from cobra import Model, Reaction, Metabolite
import os
from os.path import join
import logging


'''
Inputs:
    cfile: (str) file_path to a CSV file with compounds
Output:
    file_list: (list) dimension 1 list with strings for each line in the compounds file.
'''
def compounds_parser(cfile):
    f = open(cfile, "r")
    file_str = f.read()
    f.close()
    file_list = file_str.split('\n')
    file_list = file_list[1:]
    return file_list


'''
Inputs:
    rfile: (str) file_path to a CSV file with reactions
Output:
    file_list: (list) dimension 1 list with strings for each line in the reactions file.
'''
def reactions_parser(rfile):
    f = open(rfile, "r")
    file_str = f.read()
    f.close()
    file_list = file_str.split('\n')
    file_list = file_list[2:]
    return file_list


'''
Inputs:
    c_f_list: (list) A list made from the compounds file, each element is a row from the file

Outputs:
    compound_dict: (dict) A dictionary mapping compound names to metabolites in cobrapy format.

'''
def add_compounds_to_compound_dict(c_f_list):

    compound_dict = dict()

    for line in c_f_list:
        compound = line.split(',')
        met = Metabolite(
                compound[0], #compound id
                formula = compound[2], # formula
                name = compound[1],
                compartment = 'c'
                )

        #We add two different formats for the compounds:
        compound_dict[compound[0] + '[c0]'] = met
        compound_dict[compound[0].replace('_c0','[c0]')] = met

    logging.debug(compound_dict)

    return compound_dict


'''
Inputs:
    r_f_list: (list) A list made from the reactions file, each element is a row from the CSV file 
    compound_dict: (dict) A dictionary mapping compound names to metabolites in cobrapy format.
    c_model: (class) The cobrapy model


Outputs:
    c_model: (class) The cobrapy model

'''
def add_reactions_to_model(r_f_list, compound_dict, c_model):

    for line in r_f_list:
        reaction_list = line.split(',')
        equation = reaction_list[8]
        logging.debug(equation)
        equation_halves = split_reaction(equation)
        part_1 = turn_reaction_half_into_list(equation_halves["halves_list"][0])
        part_2 = turn_reaction_half_into_list(equation_halves["halves_list"][1])
        compound_coefficient_list = convert_compound_lists_to_proper_form(part_1, part_2, equation_halves["direction"])
        logging.debug(compound_coefficient_list)
        reaction = Reaction(reaction_list[0])
        reaction.name = reaction_list[4]
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
        for cmpnd in compound_coefficient_list:
            logging.debug(cmpnd)
            if cmpnd[1] in compound_dict:
                reaction.add_metabolites({compound_dict[cmpnd[1]] : cmpnd[0]})
            else:
                logging.debug("Did not find compound: " + cmpnd[1])
                raise Exception("Reaction incomplete, cannot continue.")
        logging.debug(reaction.reaction)
        c_model.add_reactions([reaction])

    return c_model



def turn_reaction_half_into_list(half_rxn_str):
        output_list = []
        compounds = half_rxn_str.split('+')
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



#returns a list with 2 parts.
def split_reaction(reaction_str):

    if ">" in reaction_str and "<" in reaction_str:
        logging.debug("Double reaction: " + reaction_str)
        return {'halves_list' : reaction_str.split('<=>'), "direction": "<>"}
    elif ">" in reaction_str:
        return { "halves_list": reaction_str.split("=>"), "direction" : ">"}
    elif "<" in reaction_str:
        return { "halves_list": reaction_str.split("<="), "direction" : "<"}
    else:
        raise Exception("Could not parse reaction equation into two parts")

def convert_compound_lists_to_proper_form(part_1, part_2, direction):

    if direction == ">":
        for compound_list in part_1:
            compound_list[0] = -1 * compound_list[0]
    elif direction == "<":
        for compound_list in part_2:
            compound_list[0] = -1 * compound_list[0]
    else:
        raise Exception("Cannot recognize direction.")
    return part_1 + part_2





def main():
    logging.basicConfig(level=logging.DEBUG)

    Main_Model = cobra.io.read_sbml_model('/Users/omreeg/Programs/Arkin_Lab_Research_Home/Current_Projects/Mika/FBA_Learn_Python/Examples/iAF1260.xml')
    logging.debug(len(Main_Model.metabolites))
    logging.debug(type(Main_Model.metabolites))
    x = Main_Model.metabolites.get_by_id('cpd00011')
    logging.debug(x)
    for i in range(10):
        logging.debug(Main_Model.metabolites[i])

    ''' 


    compounds_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/Current_Projects/Mika/FBA_Learn_Python/Examples/Added_Compounds_EColi_Mika_iAF1260.csv'
    reactions_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/Current_Projects/Mika/FBA_Learn_Python/Examples/Added_Reactions_EColi_Mika_iAF1260.csv'
    compounds_file_list = compounds_parser(compounds_file_name)
    reaction_file_list = reactions_parser(reactions_file_name)
    compounds_dict = add_compounds_to_compound_dict(compounds_file_list)
    Main_Model = add_reactions_to_model(reaction_file_list, compounds_dict, Main_Model)
    '''

   
    




main()
