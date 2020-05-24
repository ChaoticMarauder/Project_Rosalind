from Chapter1 import hamming_distance
from Chapter1 import neighbours

def motif_enumeration(dna, k, d):
    patterns = []
    num_dna = len(dna)
    for seq in dna:
        len_seq = len(seq)
        for i in range(len_seq-k+1):
            pattern = seq[i:i+k]
            
            neighbourhood = neighbours(pattern, d)
            
            for pattern_dash in neighbourhood:
                count=0
                for seq_comp in dna:
                    len_seq_c = len(seq_comp)
                    for j in range(len_seq_c-k+1):
                        if hamming_distance(pattern_dash, seq_comp[j:j+k])<=d:
                            count=count+1
                            break
                
                if(count == num_dna):
                    patterns.append(pattern_dash)
    
    final_patterns=[]
    for pattern in patterns:
        if pattern not in final_patterns:
            final_patterns.append(pattern)
            
    return final_patterns
                    
def main():
    with open('datasets/rosalind_ba2a.txt') as input_data:
         dna = input_data.read().strip().split('\n')
        
    k,d = dna[0].split()
    
    dna = dna[1:]
    
    k = int(k)
    d = int(d)
    

    final_patterns = motif_enumeration(dna, k, d)
    
    print(' '.join(final_patterns))
    
    with open('solutions/rosalind_ba2a.txt', 'w') as output_file:
        output_file.write(' '.join(final_patterns))
        
if(__name__=='__main__'):
    main()