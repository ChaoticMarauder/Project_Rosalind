import numpy as np
def longest_common_subseq(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    
    s = np.zeros([m+1, n+1], dtype=int)
    
    for i in range(m):
        for j in range(n):
            if seq1[i] == seq2[j]:
                s[i+1][j+1] = s[i][j] + 1
            else:
                s[i+1][j+1] = max(s[i+1][j], s[i][j+1])
        
    lcs=''
    
    i, j = m, n
    
    while(i*j!=0):
        if s[i][j] == s[i-1][j]:
            i=i-1
        elif s[i][j] == s[i][j-1]:
            j=j-1    
        else:
            if seq1[i-1] == seq2[j-1]:
                lcs = seq1[i-1] + lcs
                i=i-1
                j=j-1
                
    return lcs

def main():
    with open('datasets/rosalind_ba5c.txt') as input_file:
        seq1, seq2 = input_file.read().strip().split('\n')
            
   
    lcs = longest_common_subseq(seq1, seq2)
    
    print(lcs)
    
    
    with open('solutions/rosalind_ba5c.txt', 'w') as output_file:
        output_file.write(lcs) 
        
if __name__ == '__main__':
    main()
                