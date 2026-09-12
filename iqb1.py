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

#doing 1.c: tracing back from the max score position until we reach a score of 0
def traceback(max_row, max_col):
    alignment1 = []
    alignment2 = []
    visited_path = [(max_row, max_col)]
    current_row = max_row
    current_col = max_col
    while current_row > 0 and current_col > 0 and direction_matrix[current_row][current_col] != 0:
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

alignment1, alignment2, visited_path = traceback(max_pos[0], max_pos[1])
print("\nc.local alignment:")
print("Sequence 1:", "".join(alignment1))
print("Sequence 2:", "".join(alignment2))
print("\nvisited path from max score position:", visited_path)

#doing 1.d: Report the optimal local alignment(s), using a vertical bar (|) for a match and a hyphen (−) for a
#gap.   
print("\nd.optimal local alignment(s), using a vertical bar (|) for a match and a hyphen (−) for a gap:")
print("Sequence 1:", "".join(alignment1))
print("Sequence 2:", "".join(alignment2))

#doing 1.e: Report the optimal local alignment
print("\ne. Optimal local alignment:") 
seq1_aligned = "".join(alignment1) 
seq2_aligned = "".join(alignment2) 

match_line = "" 
for base1, base2 in zip(seq1_aligned, seq2_aligned): 
    if base1 == "-" or base2 == "-": 
        match_line += "-" 
    elif base1 == base2: 
        match_line += "|" 
    else: 
        match_line += " " 

print(seq1_aligned) 
print(match_line) 
print(seq2_aligned) 
                                         
