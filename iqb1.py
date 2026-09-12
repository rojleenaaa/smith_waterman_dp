import numpy as np
import pandas as pd

seq1 = "GATTACAAGTCC"
seq2 = "TTACAGTCA"

match = 3
mismatch_penalty = -2
gap_penalty = -2

m, n = len(seq1), len(seq2)

score_matrix = np.zeros((m+1, n+1), dtype=int)
direction_matrix = np.zeros((m+1, n+1), dtype=int)


# doing 1.a: constructing local alignment matrix

for i in range(1, m+1):

    for j in range(1, n+1):

        if seq1[i-1] == seq2[j-1]:
            diagonal_score = score_matrix[i-1][j-1] + match
        else:
            diagonal_score = score_matrix[i-1][j-1] + mismatch_penalty

        up_score = score_matrix[i-1][j] + gap_penalty
        left_score = score_matrix[i][j-1] + gap_penalty

        best_score = max(0, diagonal_score, up_score, left_score)

        score_matrix[i][j] = best_score

        if best_score == 0:
            direction_matrix[i][j] = 0

        elif best_score == diagonal_score:
            direction_matrix[i][j] = 1

        elif best_score == up_score:
            direction_matrix[i][j] = 2

        else:
            direction_matrix[i][j] = 3


score_table = pd.DataFrame(
    score_matrix,
    index=["-"] + list(seq1),
    columns=["-"] + list(seq2)
)

print("a. Local alignment matrix:")
print(score_table)


# doing 1.b: finding the maximum score and its position

max_score = np.max(score_matrix)

max_pos = np.argwhere(score_matrix == max_score)

print("\nb. Maximum score:", max_score)

print("Positions:", [tuple(pos) for pos in max_pos])


# doing 1.c: tracing back from the maximum score position
# until we reach a score of 0

def traceback(max_row, max_col):

    alignment1 = []
    alignment2 = []

    visited_path = [(max_row, max_col)]

    current_row = max_row
    current_col = max_col

    while (
        current_row > 0
        and current_col > 0
        and direction_matrix[current_row][current_col] != 0
    ):

        if direction_matrix[current_row][current_col] == 1:

            alignment1.append(seq1[current_row-1])
            alignment2.append(seq2[current_col-1])

            current_row -= 1
            current_col -= 1

        elif direction_matrix[current_row][current_col] == 2:

            alignment1.append(seq1[current_row-1])
            alignment2.append("-")

            current_row -= 1

        else:

            alignment1.append("-")
            alignment2.append(seq2[current_col-1])

            current_col -= 1

        visited_path.append((current_row, current_col))

    return alignment1[::-1], alignment2[::-1], visited_path


# Traceback from the first maximum-score position

alignment1, alignment2, visited_path = traceback(
    max_pos[0][0],
    max_pos[0][1]
)

print("\nc. Local alignment:")

print("Sequence 1:", "".join(alignment1))
print("Sequence 2:", "".join(alignment2))

print("\nVisited path from maximum score position:")
print(visited_path)


#doing 1.d: Report the optimal local alignment

print("\nd. Optimal local alignment:")

seq1_aligned = "".join(alignment1)
seq2_aligned = "".join(alignment2)
match_line = ""

for base1, base2 in zip(seq1_aligned, seq2_aligned):
    if base1 == "-" or base2 == "-":
        match_line += " "
    elif base1 == base2:
        match_line += "|"
    else:
        match_line += " "

print(seq1_aligned)
print(match_line)
print(seq2_aligned)

# doing 1.e: independently reconstructing the optimal local alignment
# using the traceback path

print("\ne. Independently reconstruct the optimal local alignment "
      "using the traceback path:")

reconstructed1 = []
reconstructed2 = []

for i in range(len(visited_path) - 1):

    row, col = visited_path[i]
    next_row, next_col = visited_path[i + 1]

    # Diagonal move

    if next_row == row - 1 and next_col == col - 1:

        reconstructed1.append(seq1[row - 1])
        reconstructed2.append(seq2[col - 1])

    # Up move

    elif next_row == row - 1:

        reconstructed1.append(seq1[row - 1])
        reconstructed2.append("-")

    # Left move

    elif next_col == col - 1:

        reconstructed1.append("-")
        reconstructed2.append(seq2[col - 1])


# Reverse because traceback moves backwards

reconstructed1.reverse()
reconstructed2.reverse()

print("Sequence 1:", "".join(reconstructed1))
print("Sequence 2:", "".join(reconstructed2))


# Independent verification of the alignment score

print("\nIndependent verification of the alignment score:")

recomputed_score = 0

for base1, base2 in zip(reconstructed1, reconstructed2):

    if base1 == "-" or base2 == "-":

        recomputed_score += gap_penalty

    elif base1 == base2:

        recomputed_score += match

    else:

        recomputed_score += mismatch_penalty


print("Alignment:", "".join(reconstructed1), "/",
      "".join(reconstructed2))

print("Recomputed score:", recomputed_score)

if recomputed_score == max_score:
    print("PASS")

else:
    print("FAIL")


# doing 1.f: why local alignment can ignore unmatched regions
# at the beginning or end of the alignment

print("\nf. Why local alignment can ignore unmatched regions "
      "at the beginning or end of the alignment?")

print(
    "\nLocal alignment is used only to find the shared regions "
    "between two sequences. Therefore, unmatched characters at "
    "the beginning or end can simply be ignored instead of being penalized."
)

print(
    "In the Smith-Waterman algorithm, the score can restart "
    "from 0, allowing the alignment to start and end whenever "
    "the best match occurs."
)