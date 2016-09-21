import numpy as np
import argparse
import random as rnd
import project1_2 as prj1


alph = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

alpha = 5
beta = 0
def gapcost(k):
    return beta + k * alpha

#         A,  C,  G,  T
mSub = [[10,  2,  5,  2],  # A
        [ 2, 10,  2,  5],  # C
        [ 5,  2, 10,  2],  # G
        [ 2,  5,  2, 10]]  # T

optimal = min

def optimal_cost(seq1, seq2, optimizer_func=max):
    global alpha, beta
    optimal = optimizer_func
    mS = [[None for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
    mD = [[None for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
    mI = [[None for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]

    # recursions
    def S(i, j):
        values = []

        if i == 0 and j == 0:
            values.append(0)
        
        if i > 0 and j > 0:
            mS_lookup = mS[i - 1][j - 1]
            char_comparison = mSub[alph[seq1[i - 1]]][alph[seq2[j - 1]]]
            if mS_lookup:
                values.append(mS_lookup + char_comparison)
            else:
                values.append(S(i - 1, j - 1) + char_comparison)

        if i > 0 and j >= 0:
            mD_lookup = mD[i][j]
            if mD_lookup:
                values.append(mD_lookup)
            else:
                values.append(D(i, j))

        if i >= 0 and j > 0:
            mI_lookup = mI[i][j]
            if mI_lookup:
                values.append(mI_lookup)
            else:
                values.append(I(i, j))

        opt_value = optimal(values)
        mS[i][j] = opt_value
        return opt_value

    def D(i, j):
        values = []

        if i > 0 and j >= 0:
            mS_lookup = mS[i - 1][j]
            if mS_lookup:
                values.append(mS_lookup - (alpha + beta))
            else:
                values.append(S(i - 1, j) - (alpha + beta))

        if i > 1 and j >= 0:
            mD_lookup = mD[i - 1][j]
            if mD_lookup:
                values.append(mD_lookup - alpha)
            else:
                values.append(D(i - 1, j) - alpha)

        opt_value = optimal(values)
        mD[i][j] = opt_value
        return opt_value

    def I(i, j):
        values = []

        if i >= 0 and j > 0:
            mS_lookup = mS[i][j - 1]
            if mS_lookup:
                values.append(mS_lookup - (alpha + beta))
            else:
                values.append(S(i, j - 1) - (alpha + beta))

        if i >=0 and j > 1:
            mI_lookup = mI[i][j - 1]
            if mI_lookup:
                values.append(mI_lookup - alpha)
            else:
                values.append(I(i, j - 1) - alpha)

        opt_value = optimal(values)
        mI[i][j] = opt_value
        return opt_value

    result = S(len(seq1), len(seq2)) 
    #print(str(backtrack(seq1, seq2, mS)) + " cost: " + str(result))
    return result

def cost(str_tuple, alpha, beta, mSub):
    def subst_cost(pair):
        cost = mSub[alph[pair[0]]][alph[pair[1]]]
        return cost

    def gap_cost(k):
        return alpha * k + beta

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

def backtrack(seq1, seq2, table):
    def subst_cost(p1, p2):
        cost = mSub[alph[p1]][alph[p2]]
        return cost

    i = len(seq1)
    j = len(seq2)
    sequence1 = ""
    sequence2 = ""

    while(i > 0 or j > 0):
        if((i > 0 and j > 0) and (table[i][j] == table[i-1][j-1] + subst_cost(seq1[i-1], seq2[j-1]))):
            i = i - 1
            j = j - 1
            sequence1 = seq1[i] + sequence1
            sequence2 = seq2[j] + sequence2

        else:
            k = 1
            while (1):
                if(i >= k and table[i][j] == table[i-k][j] - gapcost(k)):
                    
                    slashes = "-" * k
                    sequence2 = slashes + sequence2
                    sequence1 = seq1[i-k:i] + sequence1
                    i = i - k
                    break

                elif(j >= k and table[i][j] == table[i][j-k] - gapcost(k)):

                    
                    slashes = "-" * k
                    sequence1 = slashes + sequence1

                    sequence2 = seq2[j-k:j] + sequence2
                    j = j - k
                    break
                else:
                    k = k + 1
    return (sequence1, sequence2)

def generate_random_seq():
    chars = ['A', 'C', 'G', 'T']
    seq = ''
    length = rnd.randint(3,7)
    for _ in range(length):
        seq += chars[rnd.randint(0,3)]
    return seq

def bruteforce_min_cost(seq1, seq2):
    global alpha, beta, mSub
    # find all possible alignments
    alignment_tuples = prj1.optimalCostAlignment(seq1, seq2)

    # calculate cost of each, keeping track of minimum cost
    min_cost = float('inf')
    min_cost_alignment = None

    for tup in alignment_tuples:
        curr_cost = cost(tup, alpha, beta, mSub)
        print('%i   --   %s' % (curr_cost, str(tup)))
        if curr_cost < min_cost:
            min_cost = curr_cost
            min_cost_alignment = tup

    return min_cost, min_cost_alignment
    

def runTests():
    global mSub, alpha, beta, alph
    for i in range(1, 5):
    
        scoreMatrixTotal = ((alphaCost, betaCost), alphabet, scoreMatrix) = read_input_score("project_2_examples/scorematrix_1.txt")
        seq1 = read_input_fasta("project_2_examples/seq1_ex" + str(i) + ".txt")
        seq2 = read_input_fasta("project_2_examples/seq2_ex" + str(i) + ".txt")
        alpha = alphaCost
        beta = betaCost
        mSub = scoreMatrix  
        alph = alphabet

        alignment_cost = optimal_cost(seq1, seq2, min)
        print('%s\n%s\n%i\n' % (seq1, seq2, alignment_cost))

        
def testCostAgaisntBruteForce(seq1, seq2):
    allS = prj1.optimalCostAlignment(seq1, seq2)

    minCost = float("inf")

    for x in allS:
        minCost = min(minCost, cost(x, alpha, beta, mSub))

    alignment_cost = optimal_cost(seq1, seq2, min)
    print("align: " + str(alignment_cost) + " " + str(minCost))
    return alignment_cost == minCost

def project2_eval_sequences():
    alpha = -5
    beta = -5
    for i in range(1, 5):
        seq1 = read_input_fasta("project_2_examples/seq1_ex" + str(i) + ".txt")
        seq2 = read_input_fasta("project_2_examples/seq2_ex" + str(i) + ".txt")
        alignment_cost = optimal_cost(seq1, seq2, min)

    seq1 = 'tatggagagaataaaagaactgagagatctaatgtcgcagtcccgcactcgcgagatact' +\
           'cactaagaccactgtggaccatatggccataatcaaaaag'
    seq1 = seq1.upper()
    seq2 = 'atggatgtcaatccgactctacttttcctaaaaattccagcgcaaaatgccataagcacc' +\
           'acattcccttatactggagatcctccatacagccatggaa'
    seq2 = seq2.upper()
    seq3 = 'tccaaaatggaagactttgtgcgacaatgcttcaatccaatgatcgtcgagcttgcggaa' +\
           'aaggcaatgaaagaatatggggaagatccgaaaatcgaaa'
    seq3 = seq3.upper()
    seq4 = 'aaaagcaacaaaaatgaaggcaatactagtagttctgctatatacatttgcaaccgcaaa' +\
           'tgcagacacattatgtataggttatcatgcgaacaattca'
    seq4 = seq4.upper()
    seq5 = 'atgagtgacatcgaagccatggcgtctcaaggcaccaaacgatcatatgaacaaatggag' +\
           'actggtggggagcgccaggatgccacagaaatcagagcat'
    seq5 = seq5.upper()

    sequences = [seq1, seq2, seq3, seq4, seq5]
    results = [[None for _ in range(5)] for _ in range(5)]
    for i in range(5):
        for j in range(5):
            results[i][j] = optimal_cost(sequences[i], sequences[j], min)
    print(np.array(results))

    optimal_cost(seq1, seq2, min)

def testCases():
    # global mSub, alpha, beta, alph
    # main()
    # scoreMatrixTotal = ((alphaCost, betaCost), alphabet, scoreMatrix) = read_input_score('project_2_examples/scorematrix_1.txt')
    # alpha = alphaCost
    # beta = betaCost
    # mSub = scoreMatrix  
    # alph = alphabet

    assert testCostAgaisntBruteForce("ACT", "AAA")
    assert testCostAgaisntBruteForce("ACC", "ATT")
    for _ in range(3):
        assert testCostAgaisntBruteForce(generate_random_seq(), generate_random_seq())



def main():
    global mSub, alpha, beta, alph
    
    parser = argparse.ArgumentParser()

    input_seq1 = "Path to FASTA file containing the first sequence"
    parser.add_argument("--seq1", help=input_seq1)
    input_seq2 = "Path to FASTA file containing the second sequence"
    parser.add_argument("--seq2", help=input_seq2)
    input_alpha = "alpha value to be used in cost function"
    parser.add_argument("--alpha", help=input_alpha, type=int)
    input_beta = "beta value to be used in cost function"
    parser.add_argument("--beta", help=input_beta, type=int)
    input_scoreMatrix = "Path to a file containing the score matrix"
    parser.add_argument("--scoreMatrix", help=input_scoreMatrix)

    args = parser.parse_args()

    if args.seq1 and args.seq2:
        seq1 = read_input_fasta(args.seq1)
        seq2 = read_input_fasta(args.seq2)
    else:
        seq1 = read_input_fasta('fasta1.txt')
        seq2 = read_input_fasta('fasta2.txt')

    if args.scoreMatrix:
        ((alphaCost, betaCost), alphabet, scoreMatrix) = read_input_score(args.scoreMatrix)
    else:
        ((alphaCost, betaCost), alphabet, scoreMatrix) = read_input_score('project_2_examples/scorematrix_1.txt')

    # overrides whatever is in the scorematrix file
    if args.alpha:
        alphaCost = args.alpha
    if args.beta:
        betaCost = args.beta

    alpha = alphaCost
    beta = betaCost
    mSub = scoreMatrix  
    alph = alphabet
    
    # runTests()
    # print(bruteforce_min_cost('AA', 'AA'))
    # project2_eval_sequences()
    testCases()
    

if __name__ == '__main__':
    main()
    # testCases()
