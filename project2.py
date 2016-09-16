import numpy as np

import argparse




alph = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

#         A,  C,  G,  T
mSub = [[10,  2,  5,  2],  # A
        [ 2, 10,  2,  5],  # C
        [ 5,  2, 10,  2],  # G
        [ 2,  5,  2, 10]]  # T

gapcost = -5

optimal = min



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
    return read_input(aFile)[6:].replace(" ","").replace("\n", "")

def read_input(aFile):
    f = open(aFile, 'r', encoding='latin1')
    return f.read()

def read_input_score(aFile):
    with open(aFile) as f:
        lines = f.readlines()
        array = [i.split() for i in lines]
        cost = array[0][0]

        alphabet = {}

        for i in range(1, len(array)):
            alphabet[array[i][0]] = (i-1)

        scoreMatrix = [[array[i][x] for x in range(1, len(array[i]))] for i in range(1, len(array))]

        return (cost, alphabet, scoreMatrix)

def main():
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

    
    scoreMatrixTotal = (cost, alphabet, scoreMatrix) = read_input_score(args.scoreMatrix)

    #s1 = 'GGCCTAAAGGCGCCGGTCTTTCGTACCCCAAAATCTCG-GCATTTTAAGATAAGTG-AGTGTTGCGTTACACTAGCGATCTACCGCGTCTTATACT-TAAGCG-TATGCCC-AGATCTGA-CTAATCGTGCCCCCGGATTAGACGGGCTTGATGGGAAAGAACAGCTCGTC---TGTT-TAC--GTATAAACAGAATCGCCTGGGTTCGC'
    #s2 = 'GGGCTAAAGGTTAGGGTCTTTCACACTAAAGAGTGGTGCGTATCGT-GGCTAA-TGTACCGCTTC-TGGTA-TCGTGGCTTA-CG-GCCAGAC-CTACAAGTACTAGACCTGAGAACTAATCTTGTCGAGCCTTC-CATT-GA-GGG--TAATGGGAGAGAACATCGAGTCAGAAGTTATTCTTGTTTACGTAGAATCGCCTGGGTCCGC'
    #print(cost(('AATAAT', 'AA-GG-')))
    #print(cost((s1, s2)))


    


if __name__ == '__main__':
    main()
