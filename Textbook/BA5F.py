from scoring_matrices import PAM250
from Chapter5 import indel_inserted

def local_alignment(seq1, seq2, scoring_matrix, indel_penalty):
    
    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    max_score = -10
    max_i = 0
    max_j = 0
    
        
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [0, s[i-1][j]-indel_penalty, s[i][j-1]-indel_penalty, s[i-1][j-1]+scoring_matrix[seq1[i-1], seq2[j-1]]]
            s[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(s[i][j])
            
            if s[i][j] > max_score:
                max_score = s[i][j]
                max_i = i
                max_j = j
                
    a = max_i
    b = max_j
    
    aligned_seq1 = seq1[:a]
    aligned_seq2 = seq2[:b]
    
    while(a*b!=0 and backtrack_matrix[a][b]!=0):
        if backtrack_matrix[a][b] == 1:
            aligned_seq2 = indel_inserted(aligned_seq2, b)
            a = a-1
        elif backtrack_matrix[a][b] == 2:
            aligned_seq1 = indel_inserted(aligned_seq1, a)
            b = b-1
        else:
            a = a-1
            b = b-1
            
    aligned_seq1 =  aligned_seq1[a:]
    aligned_seq2 =  aligned_seq2[b:]
    
    return max_score, aligned_seq1, aligned_seq2

def main():
    
    with open('datasets/rosalind_ba5f.txt') as input_data:
        seq1, seq2 = input_data.read().strip().split('\n')

    indel_penalty = 5
    scoring_matrix = PAM250()
    max_score, aligned_seq1, aligned_seq2 = local_alignment(seq1, seq2, scoring_matrix, indel_penalty)

    print(str(max_score))
    print(aligned_seq1)
    print(aligned_seq2)
    
    with open('solutions/rosalind_ba5f', 'w') as output_data:
        output_data.write(str(max_score)+'\n')
        output_data.write(aligned_seq1+'\n')
        output_data.write(aligned_seq2)
        
if __name__ == '__main__':
    main()