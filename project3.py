import numpy as np
import argparse
import project2 as prj2
import global_linear as glin
import msa_sp_score_3k as ms_sp
import os


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

    #print(rec_helper(len(seq1), len(seq2), len(seq3)))

def backtrack(T, seq1, seq2, seq3, sp):
    i = len(seq1)
    j = len(seq2)
    k = len(seq3)

    alignment = ['', '', '']

    while i > 0 or j > 0 or k > 0:

        if i > 0 and j > 0 and k > 0 and \
        T[i][j][k] == T[i-1][j-1][k-1] + sp(seq1[i-1], seq2[j-1], seq3[k-1]):
            alignment[0] = seq1[i-1] + alignment[0]
            alignment[1] = seq2[j-1] + alignment[1]
            alignment[2] = seq3[k-1] + alignment[2]
            i -= 1
            j -= 1
            k -= 1
        elif i > 0 and j > 0 and k >= 0 and \
        T[i][j][k] == T[i - 1][j - 1][k] + sp(seq1[i - 1], seq2[j - 1], '-'):
            alignment[0] = seq1[i-1] + alignment[0]
            alignment[1] = seq2[j-1] + alignment[1]
            alignment[2] = '-' + alignment[2]
            i -= 1
            j -= 1
        elif i > 0 and j >= 0 and k > 0 and \
        T[i][j][k] == T[i - 1][j][k - 1] + sp(seq1[i - 1], '-', seq3[k - 1]):
            alignment[0] = seq1[i-1] + alignment[0]
            alignment[1] = '-' + alignment[1]
            alignment[2] = seq3[k-1] + alignment[2]
            i -= 1
            k -= 1
        elif i >= 0 and j > 0 and k > 0 and \
        T[i][j][k] == T[i][j - 1][k - 1] + sp('-', seq2[j - 1], seq3[k - 1]):
            alignment[0] = '-' + alignment[0]
            alignment[1] = seq2[j-1] + alignment[1]
            alignment[2] = seq3[k-1] + alignment[2]
            j -= 1
            k -= 1
        elif i > 0 and j >= 0 and k >= 0 and \
        T[i][j][k] == T[i - 1][j][k] + sp(seq1[i - 1], '-', '-'):
            alignment[0] = seq1[i-1] + alignment[0]
            alignment[1] = '-' + alignment[1]
            alignment[2] = '-' + alignment[2]
            i -= 1
        elif i >= 0 and j > 0 and k >= 0 and \
        T[i][j][k] == T[i][j - 1][k] + sp('-', seq2[j - 1], '-'):
            alignment[0] = '-' + alignment[0]
            alignment[1] = seq2[j-1] + alignment[1]
            alignment[2] = '-' + alignment[2]
            j -= 1
        elif i >= 0 and j >= 0 and k > 0 and \
        T[i][j][k] == T[i][j][k - 1] + sp('-', '-', seq3[k - 1]):
            alignment[0] = '-' + alignment[0]
            alignment[1] = '-' + alignment[1]
            alignment[2] = seq3[k-1] + alignment[2]
            k -= 1
        else:
            print("ERROR IN BACKTRACKING")

    print(alignment)
    assert(len(alignment[0]) == len(alignment[1]) == len(alignment[2]))
    alignment_score = 0
    for i in range(len(alignment[0])):
        alignment_score += sp(alignment[0][i], alignment[1][i], alignment[2][i])
    print('alignment score: %i' % alignment_score)


def load_brca1_globalalignments():
    sequences = [None for _ in range(8)]
    alignments = {}
    path = 'pairwise-alignments-brca1full/'
    for filename in os.listdir(path):
        f = open(path + filename, 'r')
        score = int(f.readline()[2:])
        seqi_nr = int(f.readline()[4:])
        seqi = f.readline()
        f.readline()  # skip blank line
        seqj_nr = int(f.readline()[4:])
        seqj = f.readline()
        f.readline()  # skip blank line
        f.readline()  # skip '>seqi - backtracked alignment' line
        seqi_alignment = f.readline()
        f.readline()  # skip blank line
        f.readline()  # skip '>seqj - backtracked alignment' line
        seqj_alignment = f.readline()

        if sequences[seqi_nr - 1] is None:
            sequences[seqi_nr - 1] = seqi
        if sequences[seqj_nr - 1] is None:
            sequences[seqj_nr - 1] = seqj
        alignments[(seqi_nr, seqj_nr)] = (seqi_alignment, seqj_alignment, score)

    return sequences, alignments


def sp_exact_3(seq1, seq2, seq3, substmatrix, gapcost, alphabet, incl_backtrack=False):
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
    if incl_backtrack:
        backtrack(T, seq1, seq2, seq3, sp)
    return T[last1][last2][last3]

def global_align_all_combinations(sequences_names_tuples, a, sm, g):
    sequences = [e[1] for e in sequences_names_tuples]
    sequences_names = [e[0] for e in sequences_names_tuples]

    for i, seqi in enumerate(sequences):
        for j, seqj in enumerate(sequences):
            if j <= i:
                continue  # we only want one alignment per unique pair i,j
            seqij_cost, seqij_alignment = glin.optimal_cost(seqi, seqj, a, sm, g, True)
            f = open('pairwise-alignments-brca1full/seq%i-seq%i.txt' % (i+1, j+1), 'w')
            content = '; %i\n' % seqij_cost
            content += '>seq%i\n' % (i+1)
            content += '%s\n\n' % seqi
            content += '>seq%i\n' % (j+1)
            content += '%s\n\n' % seqj
            content += '>seq%i - backtracked alignment\n' % (i+1)
            content += '%s\n\n' % seqij_alignment[0]
            content += '>seq%i - backtracked alignment\n' % (j+1)
            content += '%s\n\n' % seqij_alignment[1]
            f.write(content)
            f.close()
            print('finished seq%i and seq%i' % (i+1, j+1))


def sp_approx_2(sequences_names_tuples, alphabet, substmatrix, gapcost, outputName="output.txt"):
    a = alphabet
    sm = substmatrix
    g = gapcost
    sequences = [e[1] for e in sequences_names_tuples]
    sequences_names = [e[0] for e in sequences_names_tuples]

    min_sumcost = float('inf')
    seq_min_sumcost = ''
    for i, seqi in enumerate(sequences):
        seqi_sumcost = 0
        for j, seqj in enumerate(sequences):
            if seqi == seqj:
                continue  # same sequence, skip comparison
            seqij_cost = glin.optimal_cost(seqi, seqj, a, sm, g)[0]
            seqi_sumcost += seqij_cost
            # print('sequences %i and %i had cost: %i' % (i, j, seqij_cost))
        if seqi_sumcost < min_sumcost:
            # found a sequence with a smaller sumcost of global alignments
            min_sumcost = seqi_sumcost
            seq_min_sumcost = seqi

    #print('%s with sumcost: %i' % (seq_min_sumcost, min_sumcost))
    # above takes O(k^2 * n^2)
    # k^2 for the nested for-loops through 'sequences', n^2 for global alignment
    
    Sc = seq_min_sumcost  # center string in star tree
    # Sc = sequences[0]
    M = [Sc]
    sequences.remove(Sc)  # remaining k - 1 strings

    #print('\nalignments: (Sc, Si) where Sc is center string and Si is the ith remaining string')
    for seq in sequences:
        cost, alignment = glin.optimal_cost(Sc, seq, a, sm, g, True)
        # print(str(alignment))

        M0 = M[0]
        S0 = list(alignment[0])
        Si = alignment[1]
        curr_pos_in_M0 = 0
        pos_of_Sc_in_M0 = []
        tmp_S0 = [x for x in S0]
        while tmp_S0:
            searched_char = tmp_S0.pop(0)
            if searched_char == '-':
                continue
            for i in range(curr_pos_in_M0, len(M0)):
                if M0[i] == searched_char:
                    pos_of_Sc_in_M0.append(i)
                    curr_pos_in_M0 = i + 1
                    break
        # print('M0\t\t %s' % (M0))
        # print('Sc\t\t %s' % Sc)
        # print('S0\t\t %s' % (str(''.join(S0))))
        # print('Si\t\t %s' % (Si))
        

        Mi = ''

        # inserting initial dashes as found in M0 before processing S0 from left to right
        if pos_of_Sc_in_M0[0] > 0:
            Mi += '-' * pos_of_Sc_in_M0[0]

        # processing S0 from left to right, creating Mi as we go
        curr_pos_in_Sc = 0
        for i in range(len(S0)):
            if S0[i] == '-':
                # found a dash in S0 caused by this alignment with Si,
                # we have to insert this dash in all other entries in M on
                # this position
                if curr_pos_in_Sc == len(Sc):
                    # dash is at the end of S0 in the alignment, 
                    # just append to all strings in M instead
                    M = [s + '-' for s in M]
                else:
                    dash_pos = pos_of_Sc_in_M0[curr_pos_in_Sc]
                    M = list(map(lambda s : s[:dash_pos] + '-' + s[dash_pos:], M))

                    # we have now mutated M0, we need to recalculate pos_of_Sc_in_M0
                    # by adding 1 to every 
                    pos_of_Sc_in_M0 = [e + 1 if e >= dash_pos else e for e in pos_of_Sc_in_M0]

                # and update Mi with the character in Si on this position, 
                # which is not a dash
                Mi += Si[i]
                
            else:
                # S0 contained a character, which means whatever is in Si at
                # this pos has to go into Mi at this pos too
                Mi += Si[i]

                # if there are gaps following the character we just processed
                # present in M0 already, introduced by previous alignment 
                # processings, we have to include these as well
                if curr_pos_in_Sc + 1 <= len(pos_of_Sc_in_M0) - 1:
                    curr = curr_pos_in_Sc
                    diff = pos_of_Sc_in_M0[curr + 1] - pos_of_Sc_in_M0[curr]
                    if diff > 1:
                        Mi += '-' * (diff - 1)

                curr_pos_in_Sc += 1


        M0 = M[0]

        # inserting end dashes as found in M0 after processing S0
        if len(Mi) < len(M0):
            diff = len(M0) - len(Mi)
            Mi += '-' * diff

        M.append(Mi)

    # print M as fasta
    result_str = ''
    f = open(outputName, 'w')
    i = 0;
    for s in M:
        temp_namestr = '>'+sequences_names[i]
        print(temp_namestr)
        result_str += temp_namestr + '\n'
        print(s)
        result_str += s + '\n'
        i += 1
    f.write(result_str)
    f.close()



def run_tests(sequences_names_tuples, substmatrix, alphabet, gapcost):
    # sp_exact_3(seq1, seq2, seq3, substmatrix, 5, alphabet)
    # recursive_sp_exact_3(seq1, seq2, seq3, substmatrix, 5, alphabet)
    # sequences = ['ACGT', 'ATTCT', 'CTCGA', 'ACGGT']
    # sequences_names_tuples = [('seq%i' % i, s) for i, s in enumerate(sequences)]
    # sp_approx_2(sequences_names_tuples, alphabet, substmatrix, 5)
    load_brca1_globalalignments()


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



def generated_tests(alphabet, substmatrix):

    for i in range(0, 50):
        seq1 = prj2.generate_random_seq_with_length(30)
        seq2 = prj2.generate_random_seq_with_length(30)
        seq3 = prj2.generate_random_seq_with_length(30)
        sequences_tuple = [("none", seq1), ("none", seq2), ("none", seq3)]

        sp_exact_score = sp_exact_3(seq1, seq2, seq3, substmatrix, 5, alphabet)
        sp_approx_2(sequences_tuple, alphabet, substmatrix, 5, "output/generated.fasta")
        sp_approx_score = ms_sp.compute_sp_score("output/generated.fasta")
        calculatedMax = (2*(10-1) / 10)*sp_exact_score
        assert(sp_exact_score <= sp_approx_score and sp_approx_score <= calculatedMax)
        print('%i < %i < %i' % (sp_exact_score, sp_approx_score, calculatedMax))
        print()

def read_input_fasta(aFile):
    array = []
    inputFile = prj2.read_input(aFile)

    lines = inputFile.split(">")
    lines.pop(0)

    for i in lines:
        name = i.split('\n', 1)[0]
        value = i.split('\n', 1)[1].replace("\n", "")
        array.append((name, value))

    return array


def experiment(alphabet, substmatrix):

    for i in range(10, 210, 10):
        name = "testseqs_"+ str(i) +"_3.fasta"
        sequences_names_tuples = read_input_fasta("testseqs/" + name)
        sequences = [e[1] for e in sequences_names_tuples]
        sp_exact_score = sp_exact_3(sequences[0], sequences[1], sequences[2], substmatrix, 5, alphabet)
        sp_approx_2(sequences_names_tuples, alphabet, substmatrix, 5, "output/" + name)
        sp_approx_score = ms_sp.compute_sp_score("output/"+ name)
        calculated = (2*(i-1) / i)*sp_exact_score
        print("sp_exact_3: %i, sp_approx_score: %i, should be: %i" % (sp_exact_score, sp_approx_score, calculated))

def main():
    parser = argparse.ArgumentParser()

    help_seqs = "Path to FASTA file containing the sequences"
    parser.add_argument('sequences', help=help_seqs)
    help_substmatrix = "Path to a file containing the score matrix"
    parser.add_argument("substmatrix", help=help_substmatrix)

    args = parser.parse_args()

    sequences_names_tuples = read_input_fasta(args.sequences)
    alphabet, substmatrix = prj2.read_input_score(args.substmatrix)

    run_tests(sequences_names_tuples, substmatrix, alphabet, 5)



if __name__ == '__main__':
    main()
