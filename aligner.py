import numpy as np

MATCH = 10
MISMATCH = -1
GAP = -1

STOP = 0
LEFT = 1
UP = 2
DIAGONAL = 3


def smith_waterman(seq1, seq2):
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape=(row, col), dtype=  int)
    tracing_matrix = np.zeros(shape=(row, col), dtype= int)

    max_score = -1
    max_index = (-1, -1)

    for i in range(1, row):
        for j in range(1, col):
            match_value = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            diagonal_score = matrix[i - 1, j - 1] + match_value
            vertical_score = matrix[i - 1, j] + GAP
            horizontal_score = matrix[i, j - 1] + GAP

            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)

            if matrix[i, j] == 0:
                tracing_matrix[i, j] = STOP

            elif matrix[i, j] == horizontal_score:
                tracing_matrix[i, j] = LEFT

            elif matrix[i, j] == vertical_score:
                tracing_matrix[i, j] = UP

            elif matrix[i, j] == diagonal_score:
                tracing_matrix[i, j] = DIAGONAL

            if matrix[i, j] >= max_score:
                max_index = (i, j)
                max_score = matrix[i, j]

    aligned_seq1 = ""
    aligned_seq2 = ""
    current_aligned_seq1 = ""
    current_aligned_seq2 = ""
    (max_i, max_j) = max_index

    while tracing_matrix[max_i, max_j] != STOP:
        if tracing_matrix[max_i, max_j] == DIAGONAL:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1

        elif tracing_matrix[max_i, max_j] == UP:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'
            max_i = max_i - 1

        elif tracing_matrix[max_i, max_j] == LEFT:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]
            max_j = max_j - 1

        aligned_seq1 += current_aligned_seq1
        aligned_seq2 += current_aligned_seq2

    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return aligned_seq1, aligned_seq2

def needleman_wunsch(seq1, seq2, gap_penalty=-1, match_reward=3, mismatch_penalty=-1):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    dp = [[0 for j in range(len_seq2 + 1)] for i in range(len_seq1 + 1)]

    for i in range(1, len_seq1 + 1):
        dp[i][0] = i * gap_penalty
    for j in range(1, len_seq2 + 1):
        dp[0][j] = j * gap_penalty

    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            match = dp[i - 1][j - 1] + (match_reward if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(match, delete, insert)

    align1 = ""
    align2 = ""
    i, j = len_seq1, len_seq2
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + (match_reward if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            j -= 1
    return (align2)
