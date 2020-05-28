from rosalind import parse_fasta

def indel_inserted(seq, i):
    return seq[:i] + '-' + seq[i:]
def edit_distance_alignment(seq1, seq2):
    
    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    for i in range(1, m+1):
        s[i][0] = i
        
    for i in range(1, n+1):
        s[0][i] = i

    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [s[i-1][j]+1, s[i][j-1]+1, s[i-1][j-1] + [1,0][seq1[i-1]==seq2[j-1]]]
            s[i][j] = min(score)
            backtrack_matrix[i][j] = score.index(s[i][j])
            
    i = m
    j = n
    
    aligned_seq1 = seq1
    aligned_seq2 = seq2
    
    min_score = s[i][j]
    
    while(i*j!=0):
        if backtrack_matrix[i][j] == 0:
            i = i-1
            aligned_seq2 = indel_inserted(aligned_seq2, j)
        
        elif backtrack_matrix[i][j] == 1:
            j = j-1
            aligned_seq1 = indel_inserted(aligned_seq1, i)
            
        else:
            i = i-1
            j = j-1
            
    for indel in range(i):
        aligned_seq1 = indel_inserted(aligned_seq1, indel)
    for indel in range(j):
        aligned_seq2 = indel_inserted(aligned_seq2, indel)
        
    return min_score, aligned_seq1, aligned_seq2

def main():
    
    file_name='datasets/rosalind_edta.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
    
    seq1 = seq_list[0]
    seq2 = seq_list[1]

   
    min_score, aligned_seq1, aligned_seq2 = edit_distance_alignment(seq1, seq2)

    print(str(min_score))
    print(aligned_seq1)
    print(aligned_seq2)
    
    with open('solutions/rosalind_edta', 'w') as output_data:
        output_data.write(str(min_score)+'\n')
        output_data.write(aligned_seq1+'\n')
        output_data.write(aligned_seq2)
        
if __name__ == '__main__':
    main()