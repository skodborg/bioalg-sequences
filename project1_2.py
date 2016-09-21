
import numpy as np

#seq1 = "AATAAT"
#seq2 = "AAGG"

#seq1 = "TCCAGAGA"
#seq2 = "TCGAT"

seq1 = "CGTGTCAAGTCT"
seq2 = "ACGTCGTAGCTAGG"

#gapCost = -5

'''mSub = [[10, 2, 5, 2],
		 [2, 10, 2, 5],
		 [5, 2, 10, 2],
		 [2, 5, 2, 10]]'''

gapCost = 0
mSub = [[0, 0, 0, 0],
		[0, 0, 0, 0],
		[0, 0, 0, 0],
		[0, 0, 0, 0]]


d = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}


def optimalCostAlignment(sequence1, sequence2):
	table = [[0 for x in range(len(sequence2) + 1)] for x in range(len(sequence1) + 1)]

	for x in range(1, len(table)):
		table[x][0] = gapCost * x

	for x in range(1, len(table[0])):
		table[0][x] = gapCost * x



	for y in range(1, len(table)):
		for x in range(1, len(table[0])):
			above = table[y-1][x] + gapCost
			left = table[y][x-1] + gapCost

			char1 = sequence1[y-1]
			char2 = sequence2[x-1]

			diagonal = table[y-1][x-1] + mSub[d[char1]][d[char2]]

			cost = max(above, left, diagonal)

			table[y][x] = cost

	table = np.array(table)


	backtrackingSequences = []

	def backTracking(sequences, x, y):
		result = table[y][x]
		if(x == 0 and y == 0):
			backtrackingSequences.append(sequences)
			return


		if(y > 0 and table[y-1][x] + gapCost == result):
			tempSeq1 = sequence1[y-1] + sequences[0]
			tempSeq2 = "-" + sequences[1]
			backTracking((tempSeq1, tempSeq2), x, y-1)


		if(x > 0 and table[y][x-1] + gapCost == result):
			tempSeq1 = "-" + sequences[0]
			tempSeq2 = sequence2[x-1] + sequences[1]
			backTracking((tempSeq1, tempSeq2), x-1, y)


		char1 = sequence1[y-1]
		char2 = sequence2[x-1]

		if((y > 0 and x > 0) and table[y-1][x-1] + mSub[d[char1]][d[char2]] == result):
			tempSeq1 = char1 + sequences[0]
			tempSeq2 = char2 + sequences[1]

			backTracking((tempSeq1, tempSeq2), x-1, y-1)


	backTracking(("",""), len(table[0])-1, len(table)-1)
	return backtrackingSequences


def read_input(aFile):
    f = open(aFile, 'r', encoding='latin1')
    return f.read()




def main():
	#fastaSequence1 = read_input("project_2_examples/seq1_ex1.txt")[6:].replace(" ","").replace("\n", "").upper()
	#fastaSequence2 = read_input("project_2_examples/seq2_ex1.txt")[6:].replace(" ","").replace("\n", "").upper()
	backtrack = optimalCostAlignment("AA", "AA")
	print(len(backtrack))


if __name__ == '__main__':
    main()


