import numpy as np

def print_alignments(alignment_result):
    (alignment_score, alignments) = alignment_result
    print(alignment_score)
    for (alignment_a, alignment_b) in alignments:
        print(alignment_a)
        print(alignment_b)
        print()

def needleman_wunsch(sequence_a, sequence_b, match_score, mismatch_score, gap_score):
    """
    performs needleman wunsch sequence alignment as described in:
    Needleman, Saul B. & Wunsch, Christian D. (1970).
    "A general method applicable to the search for similarities in the amino acid sequence of two proteins".
    Journal of Molecular Biology.
    Args:
    sequence a, sequence b: the two strings to be aligned
    match_score: the score given to a match between two equal characters
    mismatch_score: the score given to a match between two unequal characters
    gap_score: the score given to a match between a character and a gap
    """
    # number of rows in the dynamic programming tables
    m = len(sequence_a) + 2
    # number of columns in the dynamic programming tables
    n = len(sequence_b) + 2
    # the dynamic programming tables
    score_table = np.array([[None] * n for i in range(m)])
    arrow_table = np.array([[None] * n for i in range(m)])
    # the output list of alignments
    alignments = []

    def initialize():
        """
        initializes the dynamic programming tables
        """
        # initialize the score table with sequence headers
        score_table[0, 2:] = list(sequence_b)
        score_table[2:, 0] = list(sequence_a)
        # initialize the score table with base cases
        score_table[1, 1:] = [-i for i in range(m - 1)]
        score_table[1:, 1] = [-i for i in range(n - 1)]
        # initialize the arrow table with base cases
        arrow_table[1, 2:] = [list("←") for i in range(n - 2)]
        arrow_table[2:, 1] = [list("↑") for i in range(m - 2)]
    
    def build_tables():
        """
        populates the dynamic programming tables
        """
        for i in range(2, m):
            for j in range(2, n):
                # insert a gap in sequence_a
                up_score = score_table[i - 1, j] + gap_score
                # insert a gap in sequence_b
                left_score = score_table[i, j - 1] + gap_score
                if score_table[i, 0] == score_table[0, j]:
                    # match equal characters
                    diagonal_score = score_table[i - 1, j - 1] + match_score
                else:
                    # match unequal characters
                    diagonal_score = score_table[i - 1, j - 1] + mismatch_score
                # determine the optimal score and store it in the score table
                max_score = max(up_score, left_score, diagonal_score)
                score_table[i, j] = max_score
                # determine which arrows describe the optimal step, not necesssarily unique
                arrows = []
                if up_score == max_score:
                    arrows += "↑"
                if left_score == max_score:
                    arrows += "←"
                if diagonal_score == max_score:
                    arrows += "↖"
                # put the list of arrows into the arrow table
                arrow_table[i, j] = arrows

    def find_alignments():
        """
        use the arrow_table to construct the alignments
        we start at the bottom right and follow the arrows back to (1,1)
        a stack is used to facilitate backtracking so that all optimal alignments are found
        """
        alignment_a = []
        alignment_b = []
        i = m - 1
        j = n - 1
        stack = []

        while True:
            # if we reach (1,1) then the current alignment is complete
            if (i, j) == (1, 1):
                # current alignment is complete so add it to the list of alignments
                alignments.append(("".join(alignment_a), "".join(alignment_b)))
                # check if there are any backtracking states in the stack
                if stack == []:
                    # if not then then all alignments have been found and we are done
                    break
                else:
                    # stack is not empty so we backtrack
                    (i, j, arrows, alignment_a, alignment_b) = stack.pop()
            else:
                # general case, copy arrows for this step from from arrow_table
                arrows = list(arrow_table[i, j])
            # process first arrow for this step
            # push the remaining arrows, if any, on to the stack for later backtracking
            arrow = arrows.pop()
            if arrows != []:
                # save this state for later backtracking
                stack.append((i, j, arrows, list(alignment_a), list(alignment_b)))
            # process current arrow
            if arrow == "↖":
                # match the current position in both sequences
                alignment_a.insert(0, sequence_a[i - 2])
                alignment_b.insert(0, sequence_b[j - 2])
                # follow the arrow with a diagonal step
                i = i - 1
                j = j - 1
            elif arrow == "↑":
                # match current position in sequence a with a gap in sequence b
                alignment_b.insert(0, " ")
                alignment_a.insert(0, sequence_a[i - 2])
                # follow the arrow upwards
                i = i - 1
            elif arrow == "←":
                # match current position in sequence b with a gap in sequence a
                alignment_b.insert(0, sequence_b[j - 2])
                alignment_a.insert(0, " ")
                # follow the arrow to the left
                j = j - 1

    initialize()
    build_tables()
    find_alignments()
    # the optimal alignment score is found in the bottom right of the score table
    alignment_score = score_table[m - 1, n - 1]
    return (alignment_score, alignments)

# support use via command line
if __name__ == '__main__':
    import sys
    alignment_result = needleman_wunsch(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
    print_alignments(alignment_result)