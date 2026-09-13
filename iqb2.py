import numpy as np
import pandas as pd

seq1 = "GATTACAAGTCC"
seq2 = "TTACAGTCA"

match = 2
mismatch_penalty = -3
gap_penalty = -1

m, n = len(seq1), len(seq2)

score_matrix = np.zeros((m+1, n+1), dtype=int)
#The score matrix stores How good is the alignment at this position

direction_matrix = np.zeros((m+1, n+1), dtype=int)
#the directional matrix stores where did this score come from

#doing 1.a initialising the scoring matrix 

for i in range(m+1): #iterating through ecah row
    score_matrix[i][0] = i * gap_penalty
for j in range(n+1):
    score_matrix[0][j] = j * gap_penalty

#doing 1.b: constructionh the complete global-alignment scoring matrix

for i in range(1,m+1):
    for j in range(1,n+1):
        if seq1[i-1] == seq2[j-1]:
            diagonal_score =score_matrix[i-1][j-1] + match
        else:
            diagonal_score = score_matrix[i-1][j-1] + mismatch_penalty

        up_score = score_matrix[i-1][j] + gap_penalty

        left_score = score_matrix[i][j-1] + gap_penalty

        best_score = max(diagonal_score,up_score,left_score)
        score_matrix[i][j] = best_score

        if best_score == diagonal_score:
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

print("b. global alignment matrix:")
print(score_table)

#doing 1.c: Report the optimal global alignment
curr_row, curr_col = i-1, j-1
aligned_seq1 = ""
aligned_seq2 = ""

while curr_row > 0 or curr_col >0:
    if curr_row ==0:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[curr_col - 1] + aligned_seq2
        curr_col -= 1
    elif curr_col == 0:
        # Only seq1 bases remain -> must be gaps in seq2
        aligned_seq1 = seq1[curr_row - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        curr_row -= 1
    elif direction_matrix[curr_row][curr_col] == 1:
        aligned_seq1 = seq1[curr_row - 1] + aligned_seq1
        aligned_seq2 = seq2[curr_col - 1] + aligned_seq2
        curr_row -= 1
        curr_col -= 1
    elif direction_matrix[curr_row][curr_col] == 2:
        aligned_seq1 = seq1[curr_row - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        curr_row -= 1
    else:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[curr_col - 1] + aligned_seq2
        curr_col -= 1
 
match_line = ""
for base1, base2 in zip(aligned_seq1, aligned_seq2):
    if base1 == "-" or base2 == "-":
        match_line += "-"
    elif base1 == base2:
        match_line += "|"
    else:
        match_line += " "  # mismatch
 
final_score = int(score_matrix[i - 1][j - 1])
 
print("\nc. Optimal global alignment:")
print(aligned_seq1)
print(match_line)
print(aligned_seq2)
print("Alignment score (from bottom-right cell of matrix):", final_score)
 
# Independent check, recomputed directly from the aligned strings
recomputed_score = 0
for base1, base2 in zip(aligned_seq1, aligned_seq2):
    if base1 == "-" or base2 == "-":
        recomputed_score += gap_penalty
    elif base1 == base2:
        recomputed_score += match
    else:
        recomputed_score += mismatch_penalty
print("Recomputed score from the alignment itself:", recomputed_score,
      "| matches?", recomputed_score == final_score)
 
