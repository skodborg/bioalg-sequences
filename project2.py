import numpy as np
import argparse


alph = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

alpha = 5
beta = 0
def gapcost(k, alpha=0):
    return alpha + k * beta

#         A,  C,  G,  T
mSub = [[10,  2,  5,  2],  # A
        [ 2, 10,  2,  5],  # C
        [ 5,  2, 10,  2],  # G
        [ 2,  5,  2, 10]]  # T


optimal = min

def optimal_cost(seq1, seq2):
    mS = [[None for _ in range(len(seq2))] for _ in range(len(seq1))]
    mD = [[None for _ in range(len(seq2))] for _ in range(len(seq1))]
    mI = [[None for _ in range(len(seq2))] for _ in range(len(seq1))]

    # recursions
    def S(i, j):
        # memoization, look up before calculating
        # TODO

        values = []
        # base case
        if i == 0 and j == 0:
            values.append(0)
        
        if i > 0 and j > 0:
            values.append(S(i - 1, j - 1) + mSub[alph[seq1[i - 1]]][alph[seq2[j - 1]]])

        if i > 0 and j >= 0:
            values.append(D(i, j))

        if i >= 0 and j > 0:
            values.append(I(i, j))

        return max(values)

    def D(i, j):
        values = []

        if i > 0 and j >= 0:
            values.append(S(i - 1, j) - (alpha + beta))

        if i > 1 and j >= 0:
            values.append(D(i - 1, j) - alpha)

        return max(values)

    def I(i, j):
        values = []

        if i >= 0 and j > 0:
            values.append(S(i, j - 1) - (alpha + beta))

        if i >=0 and j > 1:
            values.append(I(i, j - 1) - alpha)

        return max(values)

    result = S(len(seq1), len(seq2)) 
    return result


def cost(str_tuple):
    def subst_cost(pair):
        cost = mSub[alph[pair[0]]][alph[pair[1]]]
        return cost

    def gap_cost(k):
        return -5 * k

    str1 = str_tuple[0]
    str2 = str_tuple[1]
    costsum = 0
    curr_gapcount = 0

    for i in range(len(str1)):
        if '-' in (str1[i], str2[i]):
            # found a gap
            curr_gapcount += 1
        else:
            # found a substitution
            if curr_gapcount > 0:
                # if we have just ended a gap block, calculate cost and reset
                costsum += gap_cost(curr_gapcount)
                curr_gapcount = 0
            # add substitution cost
            costsum += subst_cost((str1[i], str2[i]))
    if curr_gapcount > 0:
        # we ended on a gap, add its cost and reset
        costsum += gap_cost(curr_gapcount)
        curr_gapcount = 0
    return costsum

def read_input_fasta(aFile):
    return read_input(aFile)[6:].replace(" ","").replace("\n", "").upper()

def read_input(aFile):
    f = open(aFile, 'r', encoding='latin1')
    return f.read()

def read_input_score(aFile):
    with open(aFile) as f:
        lines = f.readlines()
        array = [i.split() for i in lines]
        cost = (alphaCost, betaCost) = (int(array[0][0]), int(array[0][1]))

        alphabet = {}

        for i in range(1, len(array)):
            alphabet[array[i][0]] = (i-1)

        scoreMatrix = [[int(array[i][x]) for x in range(1, len(array[i]))] for i in range(1, len(array))]

        return (cost, alphabet, scoreMatrix)

def main():
    global mSub, alpha, beta
    
    parser = argparse.ArgumentParser()

    input_seq1 = "Path to FASTA file containing the first sequence"
    parser.add_argument("seq1", help=input_seq1)
    input_seq2 = "Path to FASTA file containing the second sequence"
    parser.add_argument("seq2", help=input_seq2)
    
    input_scoreMatrix = "Path to a file containing the score matrix"
    parser.add_argument("scoreMatrix", help=input_scoreMatrix)
    

    args = parser.parse_args()


    seq1 = read_input_fasta(args.seq1)
    seq2 = read_input_fasta(args.seq2)
    

    for i in range(1, 2):
    
        scoreMatrixTotal = ((alphaCost, betaCost), alphabet, scoreMatrix) = read_input_score("project_2_examples/scorematrix_1.txt")
        seq1 = read_input_fasta("project_2_examples/seq1_ex" + str(i) + ".txt")
        seq2 = read_input_fasta("project_2_examples/seq2_ex" + str(i) + ".txt")
        alpha = alphaCost
        beta = betaCost
        mSub = scoreMatrix
        print(i)
        print(optimal_cost(seq1, seq2))
    
    


if __name__ == '__main__':
    main()
