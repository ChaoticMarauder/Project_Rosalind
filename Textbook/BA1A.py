def PatternCount(dna, pattern):
    count = 0
    len_pattern = len(pattern) 
    for i in range(len(dna)-len(pattern)+1):
        if(dna[i:i+len_pattern] == pattern):
            count+=1
    return count

def main():
    with open('datasets/rosalind_ba1a.txt') as input_file:
        dna, pattern = input_file.read().strip().split('\n')
        
    count = PatternCount(dna, pattern)
    
    print(str(count))
    
    with open('solutions/rosalind_ba1a.txt', 'w') as output_file:
        output_file.write(str(count))
        
if(__name__=='__main__'):
    main()