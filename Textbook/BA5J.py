from Chapter5 import indel_inserted
from scoring_matrices import BLOSUM62

def affine_gap_penalty_global_alignment(seq1, seq2, scoring_matrix, sigma, epsilon):
    
    m = len(seq1)
    n = len(seq2)
    
    s = [[[0 for i in range(n+1)] for j in range(m+1)] for k in range(3)]
    backtrack_matrix = [[[0 for i in range(n+1)] for j in range(m+1)] for k in range(3)]
    
    for i in range(1, m+1):
        s[0][i][0] = -sigma - (i-1)*epsilon
        s[1][i][0] = -sigma - (i-1)*epsilon
        s[2][i][0] = -sigma
    
    for j in range(1, n+1):
        s[2][0][j] = -sigma - (j-1)*epsilon
        s[1][0][j] = -sigma - (j-1)*epsilon
        s[0][0][j] = -sigma
        
    for i in range(1, m+1):
        for j in range(1, n+1):
            lower_score = [s[0][i-1][j] - epsilon, s[1][i-1][j] - sigma]
            s[0][i][j] = max(lower_score)
            backtrack_matrix[0][i][j] = lower_score.index(s[0][i][j])

            middle_score = [s[0][i][j], s[1][i-1][j-1] + scoring_matrix[seq1[i-1], seq2[j-1]], s[2][i][j]]
            s[1][i][j] = max(middle_score)
            backtrack_matrix[1][i][j] = middle_score.index(s[1][i][j])
            
            upper_score = [s[2][i][j-1] - epsilon, s[1][i][j-1] - sigma]
            s[2][i][j] = max(upper_score)
            backtrack_matrix[2][i][j] = upper_score.index(s[2][i][j])
            
    aligned_seq1 = seq1 
    aligned_seq2 = seq2

    i = m
    j = n
    
    final_score = [s[0][i][j],s[1][i][j],s[2][i][j]]
    max_score = max(final_score)
    
    backtrack_score = final_score.index(max_score)
    
    while(i*j != 0):
        if backtrack_score == 0:  
            if backtrack_matrix[0][i][j] == 1:
                backtrack_score = 1
            aligned_seq2 = indel_inserted(aligned_seq2, j)
            i = i-1

        elif backtrack_score == 1:  
            if backtrack_matrix[1][i][j] == 0:
                backtrack_score = 0
            elif backtrack_matrix[1][i][j] == 2:
                backtrack_score = 2
            else:
                i = i-1
                j = j-1

        else:  
            if backtrack_matrix[2][i][j] == 1:
                backtrack_score = 1
            aligned_seq1 = indel_inserted(aligned_seq1, i)
            j = j-1
    
    for k in range(i):
        aligned_seq2 = indel_inserted(aligned_seq2, 0)
    for k in range(j):
        aligned_seq1 = indel_inserted(aligned_seq1, 0)

    return max_score, aligned_seq1, aligned_seq2

def main():
    
    with open('datasets/rosalind_ba5j.txt') as input_data:
        seq1, seq2 = input_data.read().strip().split('\n')

    sigma = 11
    epsilon = 1
    scoring_matrix = BLOSUM62()
    max_score, aligned_seq1, aligned_seq2 = affine_gap_penalty_global_alignment(seq1, seq2, scoring_matrix, sigma, epsilon)

    print(str(max_score))
    print(aligned_seq1)
    print(aligned_seq2)
    
    with open('solutions/rosalind_ba5j', 'w') as output_data:
        output_data.write(str(max_score)+'\n')
        output_data.write(aligned_seq1+'\n')
        output_data.write(aligned_seq2)
        
if __name__ == '__main__':
    main()