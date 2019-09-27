# This file goes over the steps explained at the following video:
# https://www.youtube.com/watch?v=cOUTpqlX-Xs
# For all the following, the matrix Sigma is represented by E,
#The input matrix is 'C', so we have C = UEV.


import numpy as np
from SVD.aux_1 import *


def SVD_MIT(inp_matrix):
    
    eg_values, eg_vectors = get_eigen_values_and_eigen_vectors(inp_matrix)

    
    #TEST-
    """
    print(eg_vectors)

    print("MAGNITUDES:")
    for vctr in eg_vectors:
        print(get_magnitude_of_vector(vctr))

    """
    #TEST END - 
    

    V = np.column_stack(eg_vectors)

    print("Matrix V:")
    print(V)


    E = np.diag(eig_v_to_sing_v(eg_values))
    
    print("Matrix Sigma:")
    print(E)

    #Finding 'U', we have C*V = U*E
    #U = get_U_from(inp_matrix,V,E)
    #print(U)

    #Test: Original Matrix should be UEV
    #print( np.matmul((np.matmul(U,E)),V) )

    return 0


#Step 1 multiply the matrix by its transpose
def transpose_multiplication(inp_matrix):
    
    #Get transpose of matrix
    inp_transpose = inp_matrix.transpose()

    #Multiply original matrix by its transpose
    x = np.matmul(inp_matrix, inp_transpose)

    return x

def get_U_from(inp_matrix, V, E):


    #First, we multiply inp_matrix by V
    X = np.matmul(inp_matrix, V)

    #X is also equal to UE, UE = X, so we multiply X by E(inverse), we get U

    #We get the inverse of E
    E_inverse = np.linalg.inv(E)

    #Now we multiply X by E inverse
    U = np.matmul(X, E_inverse)

    return U
