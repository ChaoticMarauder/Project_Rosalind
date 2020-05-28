from rosalind import parse_fasta
from rosalind import indel_inserted

def overlap_alignment(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    max_score = -2*(m+n)
    max_i = 0
    max_j = 0
    
        
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [s[i-1][j]-2, s[i][j-1]-2, s[i-1][j-1] + [-2,1][seq1[i-1]==seq2[j-1]]]
            s[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(s[i][j])
            
            if i == m or j == n:
                if s[i][j] > max_score:
                    max_score = s[i][j]
                    max_i = i
                    max_j = j
                
    a = max_i
    b = max_j
    
    aligned_seq1 = seq1[:a]
    aligned_seq2 = seq2[:b]
    
    while(a*b!=0):
        if backtrack_matrix[a][b] == 0:
            aligned_seq2 = indel_inserted(aligned_seq2, b)
            a = a-1
        elif backtrack_matrix[a][b] == 1:
            aligned_seq1 = indel_inserted(aligned_seq1, a)
            b = b-1
        else:
            a = a-1
            b = b-1
            
    aligned_seq1 =  aligned_seq1[a:]
    aligned_seq2 =  aligned_seq2[b:]
    
    return max_score, aligned_seq1, aligned_seq2

def main():
    
    file_name='datasets/rosalind_oap.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
    
    seq1 = seq_list[0]
    seq2 = seq_list[1]

    max_score, aligned_seq1, aligned_seq2 = overlap_alignment(seq1, seq2)

    print(str(max_score))
    print(aligned_seq1)
    print(aligned_seq2)
    
    with open('solutions/rosalind_oap.txt', 'w') as output_data:
        output_data.write(str(max_score)+'\n')
        output_data.write(aligned_seq1+'\n')
        output_data.write(aligned_seq2)
        
if __name__ == '__main__':
    main()