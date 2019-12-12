#python
import numpy as np
import pandas as pd
import cobra.test
import cobra
from cobra import Model, Reaction, Metabolite
import os
from os.path import join
import logging as loggy

'''
Inputs:
    cf_list: (list) A list made from the compounds file, each element is a row from the file

Outputs:
    met_dict: (dict) A dictionary mapping compound names to metabolites in cobrapy format.

'''
def add_met_dict(met_dict,row):
  met = Metabolite(
    row['id'], #compound id
    formula = row['formula'], # formula
    name = row['name'],
    compartment = 'c' # in cytosol
    )
  met_dict[row['id']] = met
  loggy.debug(row['id']+' has been added')
  return True
'''
Inputs:
    rf_list: (list) A list made from the reactions file, each element is a row from the CSV file 
    met_dict: (dict) A dictionary mapping compound names to metabolites in cobrapy format.
    c_model: (class) The cobrapy model
Outputs:
    c_model: (class) The cobrapy model

'''
def add_reactions_to_model(rf_list, met_dict, c_model):
    ieq = 8;
    for line in rf_list:
        rtokens = line.split(',')
        equation = rtokens[ieq]
        loggy.debug(equation)
        equation_halves = split_reaction(equation)
        part_1 = turn_reaction_half_into_list(equation_halves["halves_list"][0])
        part_2 = turn_reaction_half_into_list(equation_halves["halves_list"][1])
        compound_coefficient_list = convert_compound_lists_to_proper_form(part_1, part_2, equation_halves["direction"])
        loggy.debug(compound_coefficient_list)
        rxn_id = line[0]
        rxn_name = line[4]
        subsystem = line[2]
        lwr_bnd = 0
        up_bnd = 1000
        rxn_dict = convert_rxn_info_to_rxn_dict(compound_coefficient_list, rxn_id, rxn_name, subsystem, lwr_bnd, up_bnd)
        for cmpnd in compound_coefficient_list:
            loggy.debug(cmpnd)
            if cmpnd[1] in compound_dict:
                reaction.add_metabolites({compound_dict[cmpnd[1]] : cmpnd[0]})
            else:
                loggy.debug("Did not find compound: " + cmpnd[1])
                raise Exception("Reaction incomplete, cannot continue.")
        loggy.debug(reaction.reaction)
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
                loggy.critical('compound not found')
        half_compound_list = output_list
        return half_compound_list
#returns a list with 2 parts.
def split_reaction(reaction_str):

    if not isinstance(reaction_str, str):
        raise Exception("Input to function must be string, instead it is: " + type(reaction_str))

    if ">" in reaction_str and "<" in reaction_str:
        loggy.debug("Double reaction: " + reaction_str)
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
    ''' 
    loggy.debug("Number of metabolites: ")
    loggy.debug(len(main_model.metabolites))
    metabolites = main_model.metabolites
    for i in range(200):
        loggy.debug(metabolites[i])
    met_name = "2agpg161_c"
    if met_name in metabolites:
        ind = metabolites.index(met_name)
        x = metabolites[ind]
    loggy.debug(x)
    loggy.debug(x.formula)
    loggy.debug(x.compartment)
    loggy.debug(x.id)
    loggy.debug(x.name)
    '''
def main():
    loggy.basicConfig(level=loggy.DEBUG)
    cfile_name = './Examples/Added_Compounds_EColi_Mika_iAF1260.csv'
    rfile_name = './Examples/Added_Reactions_EColi_Mika_iAF1260.csv'
    sfile_name = './Examples/iAF1260.xml'
    try:
      cdf = pd.read_csv(cfile_name)
    except IOError:
      loggy.error('File not found: '+cfile_name)
      rdf = pd.read_csv(rfile_name)
    except IOError:
      loggy.error('File not found: '+rfile_name)
    try:
      main_model = cobra.io.read_sbml_model(sfile_name)
    except IOError:
      loggy.error('File not found: '+sfile_name)
    # check if compound is in the model
    met_list = main_model.metabolites
    mdict = {}
    for i,row in cdf.iterrows():
      if row['id'] in met_list:
        # the compound already exists in the model
        loggy.debug('metabolite already exists in the model: '+row['id'])
        continue
      else: # new compound
        loggy.debug('creating a new metabolite: '+row['id'])
        add_met_dict(mdict,row)
    for i,row in rdf.iterrows():
        # the reaction already exists in the model
        # add reaction anyhow

    
if __name__=="__main__":
    main()
