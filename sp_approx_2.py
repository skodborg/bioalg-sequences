import project3 as prj3
import project2 as prj2
import argparse


def main():
    parser = argparse.ArgumentParser()

    help_seqs = "Path to FASTA file containing the three sequences"
    parser.add_argument('sequences', help=help_seqs)
    help_substmatrix = "Path to a file containing the score matrix"
    parser.add_argument("substmatrix", help=help_substmatrix)
    help_gapcost = 'Value for the linear gapcost'
    parser.add_argument('gapcost', type=int, help=help_gapcost)

    args = parser.parse_args()

    sequences_names_tuples = prj3.read_input_fasta(args.sequences)
    alphabet, substmatrix = prj2.read_input_score(args.substmatrix)
    gapcost = args.gapcost

    prj3.sp_approx_2(sequences_names_tuples, alphabet, substmatrix, gapcost)


if __name__ == '__main__':
    main()
