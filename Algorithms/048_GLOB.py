from scoring_matrices import BLOSUM62
from rosalind import parse_fasta

def global_alignment_score(seq1, seq2, scoring_matrix, indel_penalty):

    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    
    for i in range(1, m+1):
        s[i][0] = -i*indel_penalty
    for j in range(1, n+1):
        s[0][j] = -j*indel_penalty

    
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [s[i-1][j] - indel_penalty, s[i][j-1] - indel_penalty, s[i-1][j-1] + scoring_matrix[seq1[i-1], seq2[j-1]]]
            s[i][j] = max(score)

    
    max_score = s[m][n]
        

    return max_score

def main():
    file_name='datasets/rosalind_glob.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
    
    seq1 = seq_list[0]
    seq2 = seq_list[1]
            
    indel_penalty = 5
    scoring_matrix = BLOSUM62()
    max_score = global_alignment_score(seq1, seq2, scoring_matrix, indel_penalty)
    
    print(str(max_score))
    
    
    with open('solutions/rosalind_glob.txt', 'w') as output_file:
        output_file.write(str(max_score)) 
        
if __name__ == '__main__':
    main()