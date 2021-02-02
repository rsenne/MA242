""""
the purpose of this function is to reduce an m x n matrix to reduced echelon form
#this algorithm is an implementation of the "Row Reduction Algorithm" from David C. Lay's Linear Algebra
and its Applications (5th edition pg. 15) f
"""

__author__ = "Ryan Senne"
__date__ = "Feb. 1"

import numpy as np


# this comment is left as an exercise to the programmer (I'm looking at you zachary)
def reduce_ech(mat):
    matrix = np.array(mat, dtype=np.float64)
    m, n = matrix.shape
    pivot = 0
    for r in range(m):
        if pivot >= n:
            return
        i = r
        while matrix[i][pivot] == 0:
            i += 1
            if i == m:
                i = r
                pivot += 1
                if m == pivot:
                    return
        matrix[i], matrix[r] = matrix[r], matrix[i]
        piv_v = matrix[r][pivot]
        matrix[r] = [mrx / float(piv_v) for mrx in matrix[r]]
        for i in range(m):
            if i != r:
                piv_v = matrix[i][pivot]
                matrix[i] = [iv - piv_v * rv for rv, iv in zip(matrix[r], matrix[i])]
        pivot += 1
    return matrix

# this was a test
# mtx = [
#     [1, 2, -1, -4],
#     [2, 3, -1, -11],
#     [-2, 0, -3, 22], ]
#
# b = reduce_ech(mtx)
