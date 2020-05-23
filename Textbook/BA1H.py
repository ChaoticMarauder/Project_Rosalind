from Chapter1 import hamming_distance
def approximate_pattern_matching(pattern, dna, d):
    start_positions=[]
    
    p_len=len(pattern)
    for i in range(len(dna)-p_len+1):
        pattern_dna = dna[i:i+p_len]
        if(hamming_distance(pattern, pattern_dna)<=d):
            start_positions.append(i)
            
    return start_positions

def main():
    with open('datasets/rosalind_ba1h.txt') as input_data:
        pattern, dna, d = input_data.read().strip().split('\n')
        
    d = int(d)
    start_positions=approximate_pattern_matching(pattern, dna, d)
    
    print(' '.join(list(map(str, start_positions))))
    
    with open('solutions/rosalind_ba1h.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str, start_positions))))
        
if(__name__=='__main__'):
    main()