from Chapter5 import indel_inserted

def fitting_alignment(seq1, seq2):
    
    m = len(seq1)
    n = len(seq2)
    
    f = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [f[i-1][j]-1, f[i][j-1]-1, f[i-1][j-1] + [-1,1][seq1[i-1]==seq2[j-1]]]
            f[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(f[i][j])
            
    
    i = max(enumerate([f[k][j] for k in range(n, m)]),key=lambda x: x[1])[0] + n
    j = n
    
    aligned_seq1 = seq1[:i]
    aligned_seq2 = seq2[:j]
    
    max_score = f[i][j]
    
    while(i*j != 0):
        if backtrack_matrix[i][j] == 0:
            i = i-1
            aligned_seq2 = indel_inserted(aligned_seq2, j)
        elif backtrack_matrix[i][j] == 1:
            j = j-1
            aligned_seq1 = indel_inserted(aligned_seq1, i)
        else:
            i=i-1
            j=j-1
            
    aligned_seq1 = aligned_seq1[i:]
    
    return max_score, aligned_seq1, aligned_seq2

def main():
    
    with open('datasets/rosalind_ba5h.txt') as input_data:
        seq1, seq2 = input_data.read().strip().split('\n')

   
    max_score, aligned_seq1, aligned_seq2 = fitting_alignment(seq1, seq2)

    print(str(max_score))
    print(aligned_seq1)
    print(aligned_seq2)
    
    with open('solutions/rosalind_ba5h', 'w') as output_data:
        output_data.write(str(max_score)+'\n')
        output_data.write(aligned_seq1+'\n')
        output_data.write(aligned_seq2)
        
if __name__ == '__main__':
    main()
            