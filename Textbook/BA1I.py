from Chapter1 import neighbours
from Chapter1 import NumberToPattern
from Chapter1 import PatternToNumber
from Chapter1 import approximate_pattern_matching

def frequent_words_with_mismatches(dna, k, d):
    frequent_patterns=[]
    close_array=[]
    frequency_array=[]
    
    for i in range(4**k):
        close_array.append(0)
        frequency_array.append(0)
        
    for i in range(len(dna)-k+1):
        neighbourhood = neighbours(dna[i:i+k], d)
        for pattern in neighbourhood:
            idx = PatternToNumber(pattern)
            close_array[idx] = 1
    
    for i in range(4**k):
        if close_array[i]==1:
            pattern = NumberToPattern(i,k)
            frequency_array[i] = len(approximate_pattern_matching(pattern, dna, d))

    max_value = max(frequency_array)
    for i in range(4**k):
        if frequency_array[i]==max_value:
            pattern = NumberToPattern(i,k)
            frequent_patterns.append(pattern)
            
    return frequent_patterns

def main():
    with open('datasets/rosalind_ba1i.txt') as input_data:
        dna, A = input_data.read().strip().split('\n')
        
    k, d = A.split()
    
    k = int(k)
    d = int(d)
    
    frequent_patterns = frequent_words_with_mismatches(dna, k, d)
    
    print(' '.join(frequent_patterns))
    
    with open('solutions/rosalind_ba1i.txt', 'w') as output_file:
        output_file.write(' '.join(frequent_patterns))
        
if(__name__=='__main__'):
    main()