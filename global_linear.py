import argparse
import project2 as p2

def main():
    parser = argparse.ArgumentParser()

    input_seq1 = "Path to FASTA file containing the first sequence"
    parser.add_argument("seq1", help=input_seq1)
    input_seq2 = "Path to FASTA file containing the second sequence"
    parser.add_argument("seq2", help=input_seq2)
    input_scoreMatrix = "Path to a file containing the score matrix"
    parser.add_argument("scoreMatrix", help=input_scoreMatrix)
    input_alpha = "alpha value to be used in linear cost function"
    parser.add_argument("alpha", help=input_alpha, type=int)
    input_backtrack = "returns an optimal alignment along with the cost"
    parser.add_argument('-b', '--backtrack', help=input_backtrack, action='store_true')

    args = parser.parse_args()

    seq1 = p2.read_input_fasta(args.seq1)
    seq2 = p2.read_input_fasta(args.seq2)

    alphabet, scoreMatrix = p2.read_input_score(args.scoreMatrix)

    alphaCost = args.alpha
    bool_backtrack = args.backtrack

    p2.alpha = alphaCost
    p2.beta = 0
    p2.mSub = scoreMatrix  
    p2.alph = alphabet

    cost, alignment = p2.optimal_cost(seq1, seq2, min, bool_backtrack)
    print('cost: %i' % cost)
    if bool_backtrack:
        print('>seq1')
        print(alignment[0])
        print()
        print('>seq2')
        print(alignment[1])


if __name__ == '__main__':
    main()
