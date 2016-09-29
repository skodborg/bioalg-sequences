import argparse
import project2 as prj2
import global_linear as glin


def sp_func(substmatrix, gapcost, alph):
    def hasGap(char):
        return char == '-'

    def returnFunction(char1, char2, char3):
        if hasGap(char1) and hasGap(char2) and hasGap(char3):
            return 0

        if(not hasGap(char1) and not hasGap(char2) and not hasGap(char3)):
            # return sub(A[i], B[j]) + sub(B[j], C[k]) + sub(A[i], C[k])
            return substmatrix[alph[char1]][alph[char2]] + \
                   substmatrix[alph[char2]][alph[char3]] + \
                   substmatrix[alph[char1]][alph[char3]]

        if(not hasGap(char1) and not hasGap(char2) and hasGap(char3)):
            # return sub(A[i], B[j]) + gap + gap
            return substmatrix[alph[char1]][alph[char2]] + gapcost + gapcost

        if(not hasGap(char1) and hasGap(char2) and not hasGap(char3)):
            # return gap + sub(A[i], C[k]) + gap
            return gapcost + substmatrix[alph[char1]][alph[char3]] + gapcost

        if(hasGap(char1) and not hasGap(char2) and not hasGap(char3)):
            # return gap + gap + sub(B[j], C[k])
            return gapcost + gapcost + substmatrix[alph[char2]][alph[char3]]

        return gapcost + gapcost

    return returnFunction


def recursive_sp_exact_3(seq1, seq2, seq3, substmatrix, gapcost, alphabet):
    len1 = range(len(seq1) + 1)
    len2 = range(len(seq2) + 1)
    len3 = range(len(seq3) + 1)
    sp = sp_func(substmatrix, gapcost, alphabet)
    T = [[[None for _ in len3] for _ in len2] for _ in len1]

    def rec_helper(i, j, k):
        # memoization
        if T[i][j][k]:
            return T[i][j][k]

        # base case
        if i == 0 and j == 0 and k == 0:
            return 0

        # recursive cases
        v0 = v1 = v2 = v3 = v4 = v5 = v6 = float('inf')
        if i > 0 and j > 0 and k > 0:
            rec_call = rec_helper(i - 1, j - 1, k - 1)
            v0 = rec_call + sp(seq1[i - 1], seq2[j - 1], seq3[k - 1])
        if i > 0 and j > 0 and k >= 0:
            rec_call = rec_helper(i - 1, j - 1, k)
            v1 = rec_call + sp(seq1[i - 1], seq2[j - 1], '-')
        if i > 0 and j >= 0 and k > 0:
            rec_call = rec_helper(i - 1, j, k - 1)
            v2 = rec_call + sp(seq1[i - 1], '-', seq3[k - 1])
        if i >= 0 and j > 0 and k > 0:
            rec_call = rec_helper(i, j - 1, k - 1)
            v3 = rec_call + sp('-', seq2[j - 1], seq3[k - 1])
        if i > 0 and j >= 0 and k >= 0:
            rec_call = rec_helper(i - 1, j, k)
            v4 = rec_call + sp(seq1[i - 1], '-', '-')
        if i >= 0 and j > 0 and k >= 0:
            rec_call = rec_helper(i, j - 1, k)
            v5 = rec_call + sp('-', seq2[j - 1], '-')
        if i >= 0 and j >= 0 and k > 0:
            rec_call = rec_helper(i, j, k - 1)
            v6 = rec_call + sp('-', '-', seq3[k - 1])

        min_val = min(v0, v1, v2, v3, v4, v5, v6)
        T[i][j][k] = min_val
        return min_val

    print(rec_helper(len(seq1), len(seq2), len(seq3)))


def sp_exact_3(seq1, seq2, seq3, substmatrix, gapcost, alphabet):
    len1 = range(len(seq1) + 1)
    len2 = range(len(seq2) + 1)
    len3 = range(len(seq3) + 1)
    sp = sp_func(substmatrix, gapcost, alphabet)
    T = [[[None for _ in len3] for _ in len2] for _ in len1]

    for i in len1:
        for j in len2:
            for k in len3:
                v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = float('inf')
                if i == 0 and j == 0 and k == 0:
                    v0 = 0
                if i > 0 and j > 0 and k > 0:
                    tbl_lookup = T[i - 1][j - 1][k - 1]
                    v1 = tbl_lookup + sp(seq1[i - 1], seq2[j - 1], seq3[k - 1])
                if i > 0 and j > 0 and k >= 0:
                    tbl_lookup = T[i - 1][j - 1][k]
                    v2 = tbl_lookup + sp(seq1[i - 1], seq2[j - 1], '-')
                if i > 0 and j >= 0 and k > 0:
                    tbl_lookup = T[i - 1][j][k - 1]
                    v3 = tbl_lookup + sp(seq1[i - 1], '-', seq3[k - 1])
                if i >= 0 and j > 0 and k > 0:
                    tbl_lookup = T[i][j - 1][k - 1]
                    v4 = tbl_lookup + sp('-', seq2[j - 1], seq3[k - 1])
                if i > 0 and j >= 0 and k >= 0:
                    tbl_lookup = T[i - 1][j][k]
                    v5 = tbl_lookup + sp(seq1[i - 1], '-', '-')
                if i >= 0 and j > 0 and k >= 0:
                    tbl_lookup = T[i][j - 1][k]
                    v6 = tbl_lookup + sp('-', seq2[j - 1], '-')
                if i >= 0 and j >= 0 and k > 0:
                    tbl_lookup = T[i][j][k - 1]
                    v7 = tbl_lookup + sp('-', '-', seq3[k - 1])
                T[i][j][k] = min(v0, v1, v2, v3, v4, v5, v6, v7)
    last1 = len(seq1)
    last2 = len(seq2)
    last3 = len(seq3)
    print(T[last1][last2][last3])


def run_tests(seq1, seq2, seq3, substmatrix, alphabet):
    sp_exact_3(seq1, seq2, seq3, substmatrix, 5, alphabet)
    recursive_sp_exact_3(seq1, seq2, seq3, substmatrix, 5, alphabet)


def test_sp(substmatrix, alphabet):
    #    A  C  G  T
    # A  0  5  2  5
    # C  5  0  5  2
    # G  2  5  0  5
    # T  5  2  5  0

    sp = sp_func(substmatrix, 5, alphabet)

    # Check no GAP; sub(A[i], B[j]) + sub(B[j], C[k]) + sub(A[i], C[k])
    #               sub("A", "C") + sub("C", "G") + sub("A", "G")
    #               5 + 5 + 2
    assert(sp("A", "C", "G") == 12)
    assert(sp('A', 'A', 'A') == 0)
    assert(sp('A', 'T', 'A') == 10)
    assert(sp('A', 'T', 'G') == 12)

    # Check C GAP; sub(A[i], B[j]) + gap + gap
    #              sub("G", "C") + gap + gap
    #              5 + 5 + 5
    assert(sp("G", "C", "-") == 15)
    assert(sp('A', 'A', '-') == 10)

    # Check B GAP; gap + sub(A[i], C[k]) + gap
    #              gap + sub("G", "A") + gap
    #              5 + 2 + 5
    assert(sp("G", "-", "A") == 12)
    assert(sp("G", "-", "G") == 10)

    # Check A GAP; gap + gap + sub(B[j], C[k])
    #              gap + gap + sub("C", "G")
    #              5 + 5 + 5
    assert(sp("-", "C", "G") == 15)
    assert(sp("-", "C", "C") == 10)

    # Check multiple GAP; gap + gap
    #                     5 + 5
    assert(sp("-", "C", "-") == 10)
    assert(sp("-", "-", "G") == 10)
    assert(sp("A", "-", "-") == 10)
    assert(sp("-", "-", "-") == 0)


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

    seq1 = 'GTTCCGAAAGGCTAGCGCTAGGCGCC'
    seq2 = 'ATGGATTTATCTGCTCTTCG'
    seq3 = 'TGCATGCTGAAACTTCTCAACCA'

    # seq1 = 'GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGATAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTCTCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTG'
    # seq2 = 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAACGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGA'
    # seq3 = 'CGCTGGTGCAACTCGAAGACCTATCTCCTTCCCGGGGGGGCTTCTCCGGCATTTAGGCCTCGGCGTTTGGAAGTACGGAGGTTTTTCTCGGAAGAAAGTTCACTGGAAGTGGAAGAAATGGATTTATCTGCTGTTCGAATTCAAGAAGTACAAAATGTCCTTCATGCTATGCAGAAAATCTTGGAGTGTCCAATCTGTTT'

    alphabet, substmatrix = prj2.read_input_score(args.substmatrix)

    test_sp(substmatrix, alphabet)
    run_tests(seq1, seq2, seq3, substmatrix, alphabet)


if __name__ == '__main__':
    main()
