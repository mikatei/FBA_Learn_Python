# This file serves as a place to write functions for SVD.

import numpy as np
from numpy import linalg as LA
import math

#The following is the process to get singular values for a matrix
def get_singular_values(inp_matrix):

    #The .eig function returns a tuple consisting of a vector and array. Vector contains eigenvalues.
    eigenvalues, eigenvectors = get_eigen_values_and_eigen_vectors(inp_matrix)

    #TEST:
    print(eigenvalues)


    #The singular values are the square roots of the eigenvalues
    singular_values = eig_v_to_sing_v(eigenvalues)

    return singular_values


def get_eigen_values_and_eigen_vectors(inp_matrix):

    #Get transpose of matrix
    inp_transpose = inp_matrix.transpose()

    #Multiply original matrix by its transpose
    x = np.matmul(inp_matrix, inp_transpose)

    #The .eig function returns a tuple consisting of a vector and array. Vector contains eigenvalues.
    eigenvalues, eigenvectors = LA.eig(x)

    return [eigenvalues, eigenvectors]

def get_magnitude_of_vector(inp_vctr):
        mag = 0
        for number in inp_vctr:
            mag += number*number
        return mag

def eig_v_to_sing_v(eigenvalues):

    #The singular values are the square roots of the eigenvalues
    singular_values = []

    for egv in eigenvalues:
        if egv > 0:
            singular_values.append(math.sqrt(egv))
        else:
            print("eigenvalue not greater than 0: " + str(egv))
    
    # We want the singular values to be sorted in decreasing order.
    singular_values.sort(reverse=True)

    return singular_values

