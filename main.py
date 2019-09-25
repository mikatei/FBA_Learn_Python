# This function runs FBA - this file controls all subdirectories and functions.

from stoichiometric_matrix import *
from scipy import linalg
import cobra
import sympy


def main():
    test()




def test():
    paper_file_name='/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA/Examples/paper_example.txt'
    txtbook_file_name = '/Users/omreeg/Programs/Arkin_Lab_Research_Home/FBA/Examples/textbook_94_example.txt'
    rxn_list_d2 = get_rxn_list_d2_example(txtbook_file_name)

    parsed_rxn_list_d4 = list_of_reaction_strings_to_parsed_reaction_list(rxn_list_d2)

    mtrices = create_stoichiometric_matrix(parsed_rxn_list_d4)
    S_w_cmpnds = mtrices[0]
    S = mtrices[1]

    #Creating 0 vector:
    zero_vector = []
    for i in range(len(rxn_list_d2)):
        zero_vector.append([0])
    print(zero_vector)
    b = np.array(zero_vector)

    print(S_w_cmpnds)
    print(S)
    S_r_e_form = sympy.Matrix(S).rref()
    print(S_r_e_form)

    #solving S with zero vector - problem is that S must be a square.
    #soltn = np.linalg.solve(S, b)
    #print(soltn)

main()



