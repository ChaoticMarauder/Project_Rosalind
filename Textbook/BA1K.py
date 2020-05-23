from Chapter1 import PatternToNumber

def computing_frequencies(dna, k):
    frequency_array=[]
    for i in range(4**k):
        frequency_array.append(0)
        
    for i in range(len(dna)-k+1):
        pattern=dna[i:i+k]
        index = PatternToNumber(pattern)
        frequency_array[index] = frequency_array[index]+1
        
    return frequency_array

def main():
    with open('datasets/rosalind_ba1k.txt') as input_data:
        dna, k=input_data.read().strip().split('\n')
        
    k = int(k)
    frequency_array=computing_frequencies(dna, k)
    
    print(' '.join(list(map(str, frequency_array))))
    
    with open('solutions/rosalind_ba1k.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str, frequency_array))))
        
if(__name__=='__main__'):
    main()   