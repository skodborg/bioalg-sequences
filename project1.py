import math
import numpy as np

str_alphabet = 'ACGT'
a = {}   # format: {i: C} for chars C in str_alphabet, indexed by i

#         A,  C,  G,  T
mSub = [[10,  2,  5,  2],  # A
        [ 2, 10,  2,  5],  # C
        [ 5,  2, 10,  2],  # G
        [ 2,  5,  2, 10]]  # T

gapcost = -5

optimal = max

seq1 = 'AATAAT'
seq2 = 'AAGG'

def init():
    # initialize alphabet dictionary 'a' for index lookup
    for i, c in enumerate(str_alphabet):
        a[c] = i


def main():

    init()
    print('seq1: %s\nseq2: %s' % (seq1, seq2))

    # preallocate cost table
    mCost = [[None for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]

    # base case: two empty strings takes no actions (sub, ins, del) to match
    mCost[0][0] = 0

    # fill out cost table, rows first
    for r in range(len(mCost)):
        for c in range(len(mCost[0])):
            if r == c == 0:
                # base case; already covered
                continue

            if r == 0:
                # first row; only look to the left
                mCost[r][c] = mCost[r][c - 1] + gapcost
            elif c == 0:
                # first column; only look above
                mCost[r][c] = mCost[r - 1][c] + gapcost
            else:
                # look in all three directions, decide on most optimal value
                cost_up = mCost[r - 1][c] + gapcost
                cost_left = mCost[r][c - 1] + gapcost
                cost_diag = mCost[r - 1][c - 1] + mSub[a[seq2[r - 1]]][a[seq1[c - 1]]]
                value = optimal(cost_up, cost_left, cost_diag)
                mCost[r][c] = value

    mCost_np = np.array(mCost)
    print(mCost_np)


#     ======= BACKTRACKING =======

#            A  A   T   A   A   T           
#       [  0  -5 -10 -15 -20 -25 -30]
#     A     \\
#       [ -5  10   5   0  -5 -10 -15]
#     A         \\
#       [-10   5  20==15  10   5   0]
#     G                 \\
#       [-15   0  15  22  20  15  10]
#     G                     \\
#       [-20  -5  10  17  27  25==20]


#    =========== RESULT ==========

#                AATAAT
#                AA-GG-



if __name__ == '__main__':
    main()
