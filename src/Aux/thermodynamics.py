# This file makes use of the thermodynamics equations given at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1796839/pdf/1792.pdf

import logging
import math


# The use_variable is either 1 or 0.
def test_neg_Gibbs(G_value, use_variable, K):
    if (G_value - K + K*use_variable) < 0:
        return True
    else:
        return False


#rxn_list_d3 looks like: [substrates, products, direction, reaction_name]
#substrates looks like [[num_compound, compound_name],[num_compound, compound_name], etc...
# https://en.wikipedia.org/wiki/Thermodynamic_activity
def set_Gibbs_value(rxn_list_d3):
    
    log_sum = 0

    #to calculate the log_sum, we will go through all the compounds
    substrates = rxn_list_d3[0]

    # R is the gas constant
    R = 8.31446261815324

    # T is in Kelvin - standard state temperature
    T = 298.15

    for cmpnd in substrates:
        stoich_coeff = (-1 * cmpnd[0])
        cmpnd_name = cmpnd[1]
        #how do you calculate the activity of a compound?
        
# activity is the same as concentration in our case.
def get_activity_of_compound(compound_name):

    return 0
