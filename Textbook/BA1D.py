def pattern_occurence(pattern, dna):
    list_occurences=[]
    k = len(pattern)
    for i in range(len(dna)-k):
        if(dna[i:i+k]==pattern):
            list_occurences.append(i)
            
    return list_occurences

def main():
    with open('datasets/rosalind_ba1d.txt') as input_file:
        pattern, dna = input_file.read().strip().split('\n')
        
    list_occurences = pattern_occurence(pattern, dna)
    
    print(' '.join(list(map(str, list_occurences))))
    
    with open('solutions/rosalind_ba1d.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str, list_occurences))))
        
if(__name__=='__main__'):
    main()