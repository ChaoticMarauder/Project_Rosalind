from Chapter1 import PatternToNumber
from Chapter1 import NumberToPattern

def computing_frequencies(dna, k):
    frequency_array=[]
    for i in range(4**k):
        frequency_array.append(0)
        
    for i in range(len(dna)-k+1):
        pattern=dna[i:i+k]
        index = PatternToNumber(pattern)
        frequency_array[index] = frequency_array[index]+1
        
    return frequency_array

def clump_finding(dna, k, L, t):
    clump_frequent_pattern_list=[]
    clump_array=[]
    
    for i in range(4**k-1):
        clump_array.append(0)
        
    frequency_array = computing_frequencies(dna[0:L], k)
    
    for i in range(4**k-1):
        if(frequency_array[i]>=t):
            pattern = NumberToPattern(i, k)
            clump_frequent_pattern_list.append(pattern)
    
    for i in range(1, len(dna)-L):
        first_pattern = dna[i-1:i-1+k]
        index = PatternToNumber(first_pattern)
        frequency_array[index] = frequency_array[index] - 1
        
        second_pattern = dna[i+L-k:i+L]
        index = PatternToNumber(second_pattern)
        frequency_array[index] = frequency_array[index] + 1
        
        if(frequency_array[index]>=t):
            clump_array[index]=1
            
    for i in range(4**k-1):
        if(clump_array[i]==1):
            pattern = NumberToPattern(i, k)
            clump_frequent_pattern_list.append(pattern)
            
    return clump_frequent_pattern_list


def main():
    with open('datasets/rosalind_ba1e.txt') as input_file:
        dna, A = input_file.read().strip().split('\n')
        
    k, L, t = A.split()
    
    k=int(k)
    L=int(L)
    t=int(t)
        
    clump_frequent_pattern_list = clump_finding(dna, k, L, t)
    
    print(' '.join(clump_frequent_pattern_list))
    
    with open('solutions/rosalind_ba1e.txt', 'w') as output_file:
        output_file.write(' '.join(clump_frequent_pattern_list))
        
if(__name__=='__main__'):
    main()            
