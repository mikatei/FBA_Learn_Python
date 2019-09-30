# This function runs FBA - this file controls all subdirectories and functions.

from stoichiometric_matrix import *
from scipy import linalg
import numpy as np
import cobra
import sympy
from SVD.aux_1 import get_singular_values
from Aux.aux_2 import give_upper_lower_bounds_list_d2, get_filenames
from SVD.MIT_process import SVD_MIT
from SVD.simple import quick_svd
#from linear_system import create_linear_system_from_matrix, convert_lin_sys_list_d3_to_d1_strings
from optlang_operations import stoichiomatrix_solution, model_print, make_fluxes


def main():
    test()




def test():

    filenames = get_filenames()

    print("current reactions file: " + filenames[-1])
    rxn_list_d2 = get_rxn_list_d2_example(filenames[-1])

    parsed_rxn_list_d4 = list_of_reaction_strings_to_parsed_reaction_list(rxn_list_d2)
 
    bounds_value = 5
    bounds = give_upper_lower_bounds_list_d2(parsed_rxn_list_d4, bounds_value)

    mtrices = create_stoichiometric_matrix(parsed_rxn_list_d4)
    S_w_cmpnds = mtrices[0]
    S = mtrices[1]

    print(S_w_cmpnds)
    print(S)

    
    #We assume the objective value is the index of the last added reaction.
    objective_index = len(parsed_rxn_list_d4)-1
    #The objective direction is 'max', or 'min'
    objective_direction = "max"


    #We solve the system of solutions using optlang
    status, model = stoichiomatrix_solution(S,bounds,objective_index, objective_direction)
    model_print(model)
    fluxes = np.asarray(make_fluxes(model))
    print(fluxes)

    #This is a mini-test, product_vector should be zero
    #For now product_vector should be zero
    Product_Vector = np.matmul(S,fluxes)

    print("Product Vector:")
    print(Product_Vector)
    #Here we can create the objective function



    

    '''
    lin_sys_d3 = create_linear_system_from_matrix(S)
    lin_sys_strings_d1 = convert_lin_sys_list_d3_to_d1_strings(lin_sys_d3)
    print(*lin_sys_strings_d1, sep = "\n")

    '''


main()


