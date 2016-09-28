import numpy as np
import argparse
import project2 as prj2

def sp(char1, char2, char3):
    # TODO: implement



def msa(seq1, seq2, seq3, substmatrix, gapcost):
    len1 = range(len(seq1))
    len2 = range(len(seq2))
    len3 = range(len(seq3))
    T = [[[None for _ in len3] for _ in len2] for _ in len1]

    for i in len1:
        for j in len2:
            for k in len3:
                v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float('inf')

                if i == 0 and j == 0 and k == 0:
                    v0 = 0
                if i > 0 and j > 0 and k > 0:
                    v1 = T[i - 1][j - 1][k - 1] + sp(seq1[i], seq2[j], seq3[k])
                if i > 0 and j > 0 and k >= 0:
                    v2 = T[i - 1][j - 1][k] + sp(seq1[i], seq2[j], '-')
                if i > 0 and j >= 0 and k > 0:
                    v3 = T[i - 1][j][k - 1] + sp(seq1[i], '-', seq3[k])
                if i >= 0 and j > 0 and k > 0:
                    v4 = T[i][j - 1][k - 1] + sp('-', seq2[j], seq3[k])
                if i > 0 and j >= 0 and k >= 0:
                    v5 = T[i - 1][j][k] + sp(seq1[i], '-', '-')
                if i >= 0 and j > 0 and k >= 0:
                    v6 = T[i][j - 1][k] + sp('-', seq2[j], '-')
                if i >= 0 and j >= 0 and k > 0:
                    v7 = T[i][j][k - 1] + sp('-', '-', seq3[k])
                T[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)
    
    print(np.array(T))
    
    lastx = len(seq1) - 1
    lasty = len(seq2) - 1
    lastz = len(seq3) - 1
    print(T[lastx][lasty][lastz])


def run_tests(seq1, seq2, seq3, substmatrix):
    msa(seq1, seq2, seq3, substmatrix, 5)

def main():
    parser = argparse.ArgumentParser()

    help_seq1 = "Path to FASTA file containing the first sequence"
    parser.add_argument("--seq1", help=help_seq1)
    help_seq2 = "Path to FASTA file containing the second sequence"
    parser.add_argument("--seq2", help=help_seq2)
    help_seq3 = "Path to FASTA file containing the first sequence"
    parser.add_argument("--seq3", help=help_seq3)
    help_substmatrix = "Path to a file containing the score matrix"
    parser.add_argument("substmatrix", help=help_substmatrix)

    args = parser.parse_args()

    # seq1 = read_input_fasta(args.seq1)
    # seq2 = read_input_fasta(args.seq2)
    # seq3 = read_input_fasta(args.seq3)
    seq1 = 'AAA'
    seq2 = 'CCC'
    seq3 = 'GGG'

    alphabet, substmatrix = prj2.read_input_score(args.substmatrix)

    run_tests(seq1, seq2, seq3, substmatrix)



if __name__ == '__main__':
    main()