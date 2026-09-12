import numpy as np
import pandas as pd

seq1 = "GATTACAAGTCC"
seq2 = "TTACAGTCA"

match = 3
mismatch_penalty = -2
gap_penalty = -2

m, n = len(seq1), len(seq2)
score_matrix = np.zeros((m+1, n+1), dtype = int)
direction_matrix = np.zeros((m+1, n+1), dtype = int)

#doing 1.a: constructiong local alignment matrix
for i in range(1, m+1): #iterating through each row. we start from 1 because row 0 is already initialized to 0
    for j in range(1, n+1): #iterating through each column. we start from 1 because column 0 is already initialized to 0
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

score_table = pd.DataFrame(score_matrix, index = ["-"] + list(seq1), columns = ["-"] + list(seq2))
print("a.local alignment matrix:")
print(score_table)

#doing 1.b: finding the maximum score and its position
max_score = np.max(score_matrix)
max_pos = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
print("\nb.maximum score:", max_score)
print("position:", max_pos)

