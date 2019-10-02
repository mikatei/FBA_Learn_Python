# This function runs FBA - this file controls all subdirectories and functions.

from scipy import linalg
import numpy as np
import cobra
import sympy

#The following imports require no downloads
from stoichiometric_matrix import *
from SVD.aux_1 import get_singular_values
from Aux.aux_2 import give_upper_lower_bounds_list_d2, get_filenames
from SVD.MIT_process import SVD_MIT
from SVD.simple import quick_svd
#from linear_system import create_linear_system_from_matrix, convert_lin_sys_list_d3_to_d1_strings
from optlang_operations import stoichiomatrix_solution, model_print, make_fluxes
import os





def main():
    #test()
    user_input()

    filename="myexample1.txt"
    objective_index = 0
    objective_direction = "max"
    bounds_value = 100
    PATH_TO_EXAMPLES = dir_path = os.path.dirname(os.path.realpath(__file__)) + '/Examples/'
    filepath = os.join(PATH_TO_EXAMPLES, filename)
    #quick_process(filepath,  objective_index, objective_direction, bounds_value)


def user_input():

    PATH_TO_EXAMPLES = dir_path = os.path.dirname(os.path.realpath(__file__)) + '/Examples/'
    print(PATH_TO_EXAMPLES)
    
    file_name = input("What is your example file named?  ")
    
    total_file_path = os.path.join(PATH_TO_EXAMPLES, file_name)

    is_correct = input("Is this your file location (y/n): " + total_file_path + "  " )
    while is_correct not in ["y", "n"]:
        is_correct = input("Please enter 'y' if the file location is correct, 'n' if not (lowercase)  ")
    if is_correct == "y":
        user_sub(total_file_path)
    elif is_correct == "n":
        new_filepath = input("Sorry, please input entire filepath to your example:  ")
        user_sub(new_filepath)
    else:
        print("Error, unknown")



def user_sub(total_file_path):
        bounds_string = input("Please enter an integer for the highest bound, common values are 100 or 1000:  ")
        bounds_value = int(bounds_string)
        if type(bounds_value) == int and bounds_value >= 0:
            mtrices = get_Stoichiometric_Matrix_from_File(total_file_path, bounds_value)
            S = mtrices[1]
            parsed_rxn_list_d4 = mtrices[2]
            bounds = give_upper_lower_bounds_list_d2(parsed_rxn_list_d4, bounds_value)
            n_rows = len(S[0])
            objective_index_string = input("What is the index of the reaction you want to optimize? Enter an int between 1 and " + str(n_rows) + ":  ")
            objective_index = int(objective_index_string)
            objective_direction = input("Enter 'max' if you want to maximize obj, 'min' if minimize  ")
            if type(objective_index) == int and objective_index >= 1 and objective_index <= n_rows and objective_direction in ["max", "min"]:
                   status, model = stoichiomatrix_solution(S,bounds,objective_index-1, objective_direction)
                   model_print(model)
                   fluxes = np.asarray(make_fluxes(model))
                   print(fluxes)               
                   #This is a mini-test, product_vector should be zero
                   #For now product_vector should be zero
                   Product_Vector = np.matmul(S,fluxes)
                   print("TEST: Product Vector. If all values in Product_Vector are not zero then there is an issue with the solution.")
                   print(Product_Vector)
            else:
                print("One of the inputs is incorrect. Stopping program")
                
        




def get_Stoichiometric_Matrix_from_File(filepath, bounds_value ):

    rxn_list_d2 = get_rxn_list_d2_example(filepath)

    parsed_rxn_list_d4 = list_of_reaction_strings_to_parsed_reaction_list(rxn_list_d2)
 
    bounds = give_upper_lower_bounds_list_d2(parsed_rxn_list_d4, bounds_value)

    mtrices = create_stoichiometric_matrix(parsed_rxn_list_d4)

    S_w_cmpnds = mtrices[0]
    S = mtrices[1]


    print("Stoichiometric Matrix with labelled compounds and Reactions:")
    print(S_w_cmpnds)
    print("\n")
    print("Just the Stoichiometric Matrix:")
    print(S)
    print("\n")

    return [S_w_cmpnds, S, parsed_rxn_list_d4]
    

def quick_process(filepath, objective_index, objective_direction, bounds_value):

    rxn_list_d2 = get_rxn_list_d2_example(filepath)

    parsed_rxn_list_d4 = list_of_reaction_strings_to_parsed_reaction_list(rxn_list_d2)
 
    bounds = give_upper_lower_bounds_list_d2(parsed_rxn_list_d4, bounds_value)

    mtrices = create_stoichiometric_matrix(parsed_rxn_list_d4)
    S_w_cmpnds = mtrices[0]
    S = mtrices[1]



    print("Stoichiometric Matrix with Compound Names and Reactions.")
    print(S_w_cmpnds)

    print("\n")
    print("Stoichiometric Matrix")
    print(S)
    print("\n")


    #We solve the system of solutions using optlang
    status, model = stoichiomatrix_solution(S,bounds,objective_index, objective_direction)
    model_print(model)
    fluxes = np.asarray(make_fluxes(model))

    print("Fluxes: \n")
    print(fluxes)               

    #This is a mini-test, product_vector should be zero
    #For now product_vector should be zero
    Product_Vector = np.matmul(S,fluxes)

    print("\n")
    print("Product Vector:")
    print(Product_Vector)




    return 0



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
    #objective_index = len(parsed_rxn_list_d4)-1

    #New objective indeces (subtract one from reaction number):
    objective_index = len(parsed_rxn_list_d4) - 1

    #objective_index = 10

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


    return 0






main()


