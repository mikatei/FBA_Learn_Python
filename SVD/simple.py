# Simple form of svd using numpy:

import numpy as np

def quick_svd(inp_matrix):
    u, s, vh = np.linalg.svd(inp_matrix, full_matrices=True)
    return [u,s,vh]
