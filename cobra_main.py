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
        rxn_id = line[0]
        rxn_name = line[4]
        subsystem = line[2]
        lwr_bnd = 0
        up_bnd = 1000
        rxn_dict = convert_rxn_info_to_rxn_dict(compound_coefficient_list, rxn_id, rxn_name, subsystem, lwr_bnd, up_bnd)


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


"""
Inputs:
    compound_coefficient_list: (list) of lists
Outputs:
    rxn_dict: (dict):
        rxn_id: (str)
        rxn_name: (str) Name of the reaction
        subsystem: (str) Subsystem where the reaction is meant to occur (cytosol?)
        lower_bound: (float) The lower flux bound
        upper_bound : (float) The upper flux bound
        compound_list: (list) list of compound_dicts
            compound_dict:
                compound_name: (str) Name of compound as expected to be in model_metabolites_list
                stoich_coeff: (str) int or float

"""
def convert_rxn_info_to_rxn_dict(compound_coefficient_list, rxn_id, rxn_name, subsystem, lwr_bnd,  up_bnd):
    rxn_dict = {
            "rxn_id" : rxn_id, 
            "rxn_name" : rxn_name, 
            "subsystem" : subsystem, 
            "lower_bound": lwr_bnd,
            "upper_bound": up_bnd,
            }
    compound_list = []
    for cmpnd in compound_coefficient_list:
        compound_list.append({"stoich_coeff": cmpnd[0], "compound_name": cmpnd[1]})
    rxn_dict["compound_list"] = compound_list

    return rxn_dict


'''
Inputs:
    rxn_dict: (dict):
        rxn_id: (str)
        rxn_name: (str) Name of the reaction
        subsystem: (str) Subsystem where the reaction is meant to occur (cytosol?)
        lower_bound: (float) The lower flux bound
        upper_bound : (float) The upper flux bound
        compound_list: (list) list of compound_dicts
            compound_dict:
                compound_name: (str) Name of compound as expected to be in model_metabolites_list
                stoich_coeff: (str) int or float
    model_metabolites_list: (list) A list of strings with names of compounds from the model.
        compound_object: (cobra metabolite object)
    custom_metabolites_list: (list) A list of custom metabolites:
        compound_object: (cobra metabolite object)
Output:
    main_reaction: (cobra reaction object)
'''
def make_cobra_reaction(reaction_dict, model_metabolites_list, custom_metabolites_list):

    #Preparing the metabolites:
    metabolites_to_add = []
    compound_list = reaction_dict["compound_list"]
    for cmpnd_dict in compound_list:
        compound_name = cmpnd_dict["compound_name"]
        if compound_name in custom_metabolites_list:

        

    raise Exception("Incomplete Function")
    return 0


'''
Inputs:
    cobra_model: (cobrapy model object) The Model object to test.
    compound_names: (list) A list of strings with names of metabolites expected to be in model.
Outputs:
    compounds_info_dict:
        in_model: (list) A list of strings with names of met/cmpnd in model.
        not_in_model: (list) A list of strings with names of met/cmpnd not in model.
'''
def check_if_compound_names_in_model(cobra_model, compound_names):
    metabolites_list = cobra_model.metabolites
    in_model = []
    not_in_model = []
    for name in compound_names:
        if name in metabolites_list:
            in_model.append(name)
        else:
            not_in_model.append(name)
    logging.debug("In Model: " + str(len(in_model)) + "metabolites")
    logging.debug("Not in Model: " + str(len(not_in_model)) + "metabolites")
    compounds_info_dict = {"in_model": in_model, "not_in_model": not_in_model}
    return compounds_info_dict 



"""
Inputs:
    half_rxn_str: (str)

Outputs:
    half_compound_list: (list) A list of internal lists of size 2: [[#cmpnd, nameofcmpnd],[#cmpnd,nameofcmpnd]...]

"""
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
        half_compound_list = output_list
        return half_compound_list



#returns a list with 2 parts.
def split_reaction(reaction_str):

    if not isinstance(reaction_str, str):
        raise Exception("Input to function must be string, instead it is: " + type(reaction_str))

    if ">" in reaction_str and "<" in reaction_str:
        logging.debug("Double reaction: " + reaction_str)
        return {'halves_list' : reaction_str.split('<=>'), "direction": "<>"}
    elif ">" in reaction_str:
        return { "halves_list": reaction_str.split("=>"), "direction" : ">"}
    elif "<" in reaction_str:
        return { "halves_list": reaction_str.split("<="), "direction" : "<"}
    else:
        raise Exception("Could not parse reaction equation into two parts: " + reaction_str)


"""
Inputs:
    part_1: (list) half_compound_list (a list of M*[#cmpnd, name_ofcompound])
    part_2: (list) half_compound_list (a list of N*[#cmpnd, name_ofcompound])
    direction: (str) controlled vocab: "<", ">" to represent direction of reaction.
Outputs:
    compounds_coefficients_list: (list) A list of compounds and their coefficients in the format:
        [[coeff, cmpnd_name], [coeff,cmpnd_name], ...]

"""
def convert_compound_lists_to_proper_form(part_1, part_2, direction):

    if direction == ">":
        for compound_list in part_1:
            compound_list[0] = -1 * compound_list[0]
    elif direction == "<":
        for compound_list in part_2:
            compound_list[0] = -1 * compound_list[0]
    else:
        raise Exception("Cannot recognize direction.")

    compounds_coefficients_list = part_1 + part_2
    return compounds_coefficients_list



"""
Inputs:
    filepath_to_model: (str)
Outputs:
    main_model: (cobrapy Model)
"""
def upload_sbml_model(filepath_to_model):
    main_model = cobra.io.read_sbml_model(filepath_to_model)
    return main_model


def tests():

    logging.basicConfig(level=logging.DEBUG)
   
    sbml_file = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/Current_Projects/Mika/FBA_Learn_Python/Examples/iAF1260.xml'

    main_model = upload_sbml_model(sbml_file)

    #We check if names of compounds are already in model-
    compounds_test_list = ['nadh_c','nad_c', 'co2_c', 'nh4_c','nadph_c', 'nadp_c']

    compounds_info_dict = check_if_compound_names_in_model(main_model,compounds_test_list)

    

    
    ''' 
    logging.debug("Number of metabolites: ")
    logging.debug(len(main_model.metabolites))
    metabolites = main_model.metabolites
    for i in range(200):
        logging.debug(metabolites[i])

    met_name = "2agpg161_c"
    if met_name in metabolites:
        ind = metabolites.index(met_name)
        x = metabolites[ind]
    logging.debug(x)
    logging.debug(x.formula)
    logging.debug(x.compartment)
    logging.debug(x.id)
    logging.debug(x.name)
    '''
    
    

def main():

    tests()

    ''' 


    compounds_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/Current_Projects/Mika/FBA_Learn_Python/Examples/Added_Compounds_EColi_Mika_iAF1260.csv'
    reactions_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/Current_Projects/Mika/FBA_Learn_Python/Examples/Added_Reactions_EColi_Mika_iAF1260.csv'
    compounds_file_list = compounds_parser(compounds_file_name)
    reaction_file_list = reactions_parser(reactions_file_name)
    compounds_dict = add_compounds_to_compound_dict(compounds_file_list)
    Main_Model = add_reactions_to_model(reaction_file_list, compounds_dict, Main_Model)
    '''

   
    




main()
