#Input to the calculations should be a stoichiometric matrix with size mxn.
#Output will be a list of size n in which each spot contains an equality



# S is a stoichiometric matrix, in numpy format.
def create_linear_system_from_matrix(S):

    #Each entry in the linear system will be divided into two parts [[right],[left]]
    # [right] = [[num, rxn_index],[num,rxn_index],...]. [left] will look the same.
    # right side will represent positive stoichiometric coefficients.
    lin_sys_d3 = []
    for row in S:
        eqtn_list_d2 = [[],[]]
        for i in range(len(row)):
            stoich_coeff = row[i]
            if stoich_coeff < 0:
                eqtn_list_d2[0].append([stoich_coeff, i])
            elif stoich_coeff > 0:
                eqtn_list_d2[1].append([stoich_coeff, i])
            else:
                #This means the stoichiometric coefficient is 0 - compound not involved in the 
                # reaction.
                pass
        lin_sys_d3.append(eqtn_list_d2)
    
    return lin_sys_d3


def convert_lin_sys_list_d3_to_d1_strings(lin_sys_d3):
    lin_sys_d1 = []
    for eqtn_list in lin_sys_d3:
        lin_sys_d1.append(create_linear_equation_string_from_equation(eqtn_list))
    return lin_sys_d1


# Following function takes 'eqtn_list' and turns it into a string like 3r2 + 5r3 = 4r4
#eqtn_list should look like [[right],[left]]
# [right] = [[num, index],[num,index]...]
def create_linear_equation_string_from_equation(eqtn_list):
    left_side_string = '   '
    right_side_string = '   '
    r_list = eqtn_list[0]
    #NOTE: WE ARE ADDING 1 TO THE INDEX TO MAKE IT FIT REACTION NAMES
    for coeff in r_list:
        right_side_string += str(abs(coeff[0])) + 'r' + str(coeff[1]+1) + ' + '
    l_list = eqtn_list[1]
    for coeff in l_list:
        left_side_string += str(abs(coeff[0])) + 'r' + str(coeff[1]+1) + ' + '

    
    final_string = right_side_string[:-2] + '= ' + left_side_string[:-2]
    return final_string

