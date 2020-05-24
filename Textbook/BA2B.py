import math
from Chapter1 import hamming_distance
from Chapter1 import NumberToPattern

def min_hamming_distance(pattern, dna):
    k = len(pattern)
    min_dist = len(dna)
    
    for i in range(len(dna)-k+1):
        if(hamming_distance(pattern, dna[i:i+k]) < min_dist):
            min_dist = hamming_distance(pattern, dna[i:i+k])
            
    return min_dist
        


def median_string(dna_list, k):
    distance = math.inf
    median = ''
    
    pattern_list = []
    for i in range(4**k):
        kmer = NumberToPattern(i,k)
        pattern_list.append(kmer)
        
    for pattern in pattern_list:
        score=0
        for seq in dna_list:
            score = score + min_hamming_distance(pattern, seq)
        
        if score < distance:
            distance = score
            median = pattern
            
    return median

def main():
    with open('datasets/rosalind_ba2b.txt') as input_data:
         dna = input_data.read().strip().split('\n')
        
    k = int(dna[0])
    
    dna_list = dna[1:]
    
    median = median_string(dna_list, k)
    
    print(median)
    
    with open('solutions/rosalind_ba2b.txt', 'w') as output_file:
        output_file.write(median)
        
if(__name__=='__main__'):
    main()