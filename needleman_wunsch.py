import numpy as np

def needleman_wunsch(sequence_a, sequence_b, match_score, mismatch_score, indel_score):
    m = len(sequence_a) + 2
    n = len(sequence_b) + 2
    #initialize the dynamic programming tables
    score_table = [[None] * n for i in range(m)]
    score_table[0][0] = " "
    score_table[0][1] = " "
    score_table[1][0] = " "
    for i, item in enumerate(sequence_b):
        score_table[0][i + 2] = item
    for j, item in enumerate(sequence_a):
        score_table[j + 2][0] = item
    for i in range(m - 1):
        score_table[i + 1][1] = -i
    for j in range(n - 1):
        score_table[1][j + 1] = -j
    #initialize the arrow table which describes the optimal alignments
    arrow_table = [[None] * n for i in range(m)]
    #arrow_table[1][1] = "↖"
    for i in range(1, m - 1):
        arrow_table[i + 1][1] = list("↑")
    for j in range(1, n - 1):
        arrow_table[1][j + 1] = list("←")

    # populate the dynamic programming table row by row and
    # record each optimal alignment step in the arrow table
    for i in range(2, m):
        for j in range(2, n):
            up_score = score_table[i - 1][j] + indel_score
            left_score = score_table[i][j - 1] + indel_score
            if score_table[i][0] == score_table[0][j]:
                diagonal_score = score_table[i - 1][j - 1] + match_score
            else:
                diagonal_score = score_table[i - 1][j - 1] + mismatch_score
            max_score = max(up_score, left_score, diagonal_score)
            arrows = []
            if up_score == max_score:
                arrows += "↑"
            if left_score == max_score:
                arrows += "←"
            if diagonal_score == max_score:
                arrows += "↖"
            score_table[i][j] = max_score
            arrow_table[i][j] = arrows

    #print(np.array(arrow_table))
    alignment_a = []
    alignment_b = []
    i = m - 1
    j = n - 1
    stack = []
    alignments = []
    while True:
        if (i, j) == (1, 1):
            alignments.append((alignment_a, alignment_b))
            if stack == []:
                break
            else:
                (i, j, arrows, alignment_a, alignment_b) = stack.pop()
        else:
            arrows = list(arrow_table[i][j])
        arrow = arrows.pop()
        if arrows != []:
            stack.append((i, j, arrows, list(alignment_a), list(alignment_b)))
        if arrow == "↖":
            alignment_a.insert(0, sequence_a[i - 2])
            alignment_b.insert(0, sequence_b[j - 2])
            i = i - 1
            j = j - 1
        else:
            if arrow == "↑":
                alignment_b.insert(0, " ")
                alignment_a.insert(0, sequence_a[i - 2])
                i = i - 1
            else:
                if arrow == "←":
                    alignment_b.insert(0, sequence_b[j - 2])
                    alignment_a.insert(0, " ")
                    j = j - 1
    alignment_score = score_table[m - 1][n - 1]
    return (alignment_score, alignments)

def print_alignments(alignment_result):
    (alignment_score, alignments) = alignment_result
    print(alignment_score)
    for (alignment_a, alignment_b) in alignments:
        print(alignment_a)
        print(alignment_b)
        print()
#print_alignments(needleman_wunsch(list("gattaca"), list("gcatgcu"), 1, -1, 0))
print_alignments(needleman_wunsch(list("gattaca"), list("gcatgcu"), 1, -1, 0))
#print(needleman_wunsch(list("fool"), list("foo"), 1, -1, 0))
