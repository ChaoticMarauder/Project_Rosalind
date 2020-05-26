from scoring_matrices import BLOSUM62

def global_alignment(seq1, seq2, scoring_matrix, indel_penalty):

    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]

    
    for i in range(1, m+1):
        s[i][0] = -i*indel_penalty
    for j in range(1, n+1):
        s[0][j] = -j*indel_penalty

    
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [s[i-1][j] - indel_penalty, s[i][j-1] - indel_penalty, s[i-1][j-1] + scoring_matrix[seq1[i-1], seq2[j-1]]]
            s[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(s[i][j])

    aligned_seq1 = seq1 
    aligned_seq2 = seq2

    
    max_score = s[m][n]
    i, j = m, n
    
    
    while(i*j != 0):
        if backtrack_matrix[i][j] == 0:
            aligned_seq2 = indel_inserted(aligned_seq2, j)
            i = i-1
        elif backtrack_matrix[i][j] == 1:
            aligned_seq1 = indel_inserted(aligned_seq1, i)
            j = j-1
        else:
            i = i-1
            j = j-1
            
    for indel in range(i):
        aligned_seq2 = indel_inserted(aligned_seq2, 0)
    for indel in range(j):
        aligned_seq1 = indel_inserted(aligned_seq1, 0)    

    return max_score, aligned_seq1, aligned_seq2

def indel_inserted(seq, i):
    return seq[:i] + '-' + seq[i:]


def main():
    
    with open('datasets/rosalind_ba5e.txt') as input_data:
        seq1, seq2 = input_data.read().strip().split('\n')

    indel_penalty = 5
    scoring_matrix = BLOSUM62()
    max_score, aligned_seq1, aligned_seq2 = global_alignment(seq1, seq2, scoring_matrix, indel_penalty)

    print(str(max_score))
    print(aligned_seq1)
    print(aligned_seq2)
    
    with open('solutions/rosalind_ba5e', 'w') as output_data:
        output_data.write(str(max_score)+'\n')
        output_data.write(aligned_seq1+'\n')
        output_data.write(aligned_seq2)
        
if __name__ == '__main__':
    main()