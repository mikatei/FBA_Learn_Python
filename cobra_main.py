#python
import numpy as np
import pandas as pd
import cobra.test
import cobra
from cobra import Model, Reaction, Metabolite
import os
import sys
from os.path import join
import logging as loggy
from optparse import OptionParser

MINFLUX = 0
MAXFLUX = 1000

'''
Inputs:
    met_dict: (dict) A dictionary mapping compound names to metabolites in cobrapy format.
    row: (Pandas DataFrame) A row of metabolite (compound) information
Outputs:
    met_dict: updated metabolite dictionary

'''
def add_met2dict(met_dict,row):
  new_met = Metabolite(
    row['id'], #compound id
    formula = row['formula'], # formula
    name = row['name'],
    compartment = 'c' # in cytosol
    )
  met_dict[row['id']] = new_met
  loggy.debug(row['id']+' has been added')
  return met_dict
'''
Inputs:
    met_list: (list<metabolite>) metabolites list of the model
    rxn_dict: (dict) A dictionary mapping reaction names to metabolites in cobrapy format.
    row: (Pandas DataFrame) A row of reaction information 
Outputs:
    rxn_dict: updated reaction dictionary

'''
def add_rxn2dict(met_list,rxn_dict,row):
    new_rxn = Reaction(
      row['id'], # reaction id
      name = row['name'],
      subsystem = row['kegg pathways'],
      lower_bound = MINFLUX,
      upper_bound = MAXFLUX,
    )
    rxneq = row['equation']
    eq_dict = make_eqn(met_list,rxneq)
    new_rxn.add_metabolites(eq_dict)
    rxn_dict[row['id']] = new_rxn
    loggy.debug(row['id']+' has been added')
    return rxn_dict

'''
Inputs:
    mlist: model metabolites list
    rxneq: Pandas DataFrame entry for equation
Output:
    eq_dict: add_metabolite dictionary
'''
def make_eqn(mlist,rxneq):
    eqs = rxneq.split('=')
    eq_dict = {}
    if (len(eqs) == 2):
      left_eq = eqs[0].strip('<').strip()
      eq_dict = process_half(mlist,eq_dict,left_eq,-1)
      right_eq = eqs[1].strip('>').strip()
      eq_dict = process_half(mlist,eq_dict,right_eq,1)
      return eq_dict
    else:
      loggy.error("Equation not complete: "+rxneq)
      return {}
'''
Inputs:
    mlist: model metabolites list
    eqd: (dict) add_metabolite dictionary to append the metabolites in
    half: one side of equation
    lr: -1 for LHS, 1 for RHS
Output:
    egd: updated add_metabolite dictionary

'''
def process_half(mlist,eqd,half,lr):
  for c in half.split('+'):
    cc = c.split()
    if (len(cc) == 2):
      try:
        coef = lr*float(cc[0].strip('(').strip(')').strip())
      except ValueError:
        loggy.error('Corrupted coefficient: '+ coef)
        continue
      try:
        met = mlist.get_by_id(cc[1].strip())
      except KeyError:
        loggy.error("metabolite named "+cc[1].strip()+" not found in the model")
        continue
      eqd[met] = coef
    else:
      loggy.error('Either coefficient or metabolite missing in an equation: '+ rxneqn)
      continue
  return eqd

def main(argv):
    parser = OptionParser()
    parser.add_option("-c",dest="cfile_name",help="csv file for the list of compounds")
    parser.add_option("-r",dest="rfile_name",help="csv file for the list of reactions")
    parser.add_option("-s",dest="sfile_name",help="sbml file for COBRA Model")
    parser.add_option("-d",dest="debuggy",help="True if debug mode")
    (options,args) = parser.parse_args()
    if (options.debuggy):
      loggy.basicConfig(level=loggy.DEBUG)
    '''
    cfile_name = './Examples/Added_Compounds_EColi_Mika_iAF1260.csv'
    rfile_name = './Examples/Added_Reactions_EColi_Mika_iAF1260.csv'
    sfile_name = './Examples/iAF1260.xml'
    '''
    # read all the files
    try:
      cdf = pd.read_csv(options.cfile_name)
    except IOError:
      loggy.error('File not found: '+cfile_name)
    try:
      rdf = pd.read_csv(options.rfile_name)
    except IOError:
      loggy.error('File not found: '+rfile_name)
    try:
      main_model = cobra.io.read_sbml_model(options.sfile_name)
    except IOError:
      loggy.error('File not found: '+sfile_name)
      return 0
    met_list = main_model.metabolites
    rxn_list = main_model.reactions
    mdict = {}
    rdict = {}
    for i,row in cdf.iterrows():
      if row['id'] in met_list:
        # the compound already exists in the model
        loggy.debug('metabolite already exists in the model: '+row['id'])
        continue
      else: # new compound
        loggy.debug('creating a new metabolite: '+row['id'])
        mdict = add_met2dict(mdict,row)
    # add metabolites into the model
    if (mdict):
      main_model.add_metabolites(mdict.values())
    for i,row in rdf.iterrows():
      if row['id'] in rxn_list:
        # the reaction already exists in the model
        loggy.debug('reaction already exists in the model: '+row['id'])
      else:
        add_rxn2dict(met_list,rdict,row)
    # add reaction to the model anyhow
    if (rdict):
        main_model.add_reactions(rdict.values())
    
if __name__=="__main__":
    main(sys.argv)
