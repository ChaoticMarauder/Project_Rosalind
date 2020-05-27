from rosalind import parse_fasta

def edit_distance(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    
    edit_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    for i in range(1, m+1):
        edit_matrix[i][0] = i
        
    for j in range(1, n+1):
        edit_matrix[0][j] = j
        
    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:
                edit_matrix[i][j] = edit_matrix[i-1][j-1]
            else:
                edit_matrix[i][j] = min(edit_matrix[i-1][j]+1, edit_matrix[i][j-1]+1, edit_matrix[i-1][j-1]+1)
    
    
    edit_distance_value = edit_matrix[m][n]
    
    return edit_distance_value

def main():
    file_name='datasets/rosalind_edit.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
    
    seq1 = seq_list[0]
    seq2 = seq_list[1]
            
   
    edit_distance_value = edit_distance(seq1, seq2)
    
    print(str(edit_distance_value))
    
    
    with open('solutions/rosalind_edit.txt', 'w') as output_file:
        output_file.write(str(edit_distance_value)) 
        
if __name__ == '__main__':
    main()