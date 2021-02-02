""""
this is a library of algorithms and functions relevant to linear algebra in MA242 SPR2021. I intend this to be a
learning tool for myself but, you may find this useful. Please forward any suggestions to rsenne@bu.edu
"""

__author__ = "Ryan Senne"
__date__ = "Feb. 1 2021"
__license__ = "MIT"
__version__ = "1.0.0"

import numpy as np

""""
the purpose of this function is to reduce an m x n matrix to reduced echelon form
this algorithm is an implementation of the "Row Reduction Algorithm" from David C. Lay's Linear Algebra
and its Applications (5th edition pg. 15)
"""


# this comment is left as an exercise to the programmer (I'm looking at you zachary)
def reduce_ech(mat):
    matrix = np.array(mat, dtype=np.float64)  # can swap to 128 bit for more accuracy
    m, n = matrix.shape  # get shape of input
    pivot = 0  # initialize pivot index
    # makes sure pivot index is not greater than number of columns
    # uses row swapping and self division of pivot for reduction
    for r in range(m):
        if pivot >= n:
            return
        i = r
        # if current indexed value is 0, move on
        while matrix[i][pivot] == 0:
            i += 1
            if i == m:
                i = r
                pivot += 1
                if m == pivot:
                    return
        matrix[i], matrix[r] = matrix[r], matrix[i]  # swap Rr and Ri
        piv_v = matrix[r][pivot]  # value of pivot position
        matrix[r] = [mrx / float(piv_v) for mrx in matrix[r]]  # divide pivot position by itself to get leading 1
        # perform elementary row operation for when rows are not the same
        for i in range(m):
            if i != r:
                piv_v = matrix[i][pivot]
                matrix[i] = [iv - piv_v * rv for rv, iv in zip(matrix[r], matrix[i])]  # Ri - pivot * Rr
        pivot += 1
    return matrix
