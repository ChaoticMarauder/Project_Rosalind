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