from Chapter1 import PatternCount

def frequent_words(dna, k):
    dict_kmer={}
    
    for i in range(len(dna)-k):
        kmer = dna[i:i+k]
        if kmer not in dict_kmer:
            dict_kmer[kmer] = PatternCount(dna, kmer)
            
    max_value=0
    
    for key, value in dict_kmer.items():
        if(value >= max_value):
            max_value = value
            
    list_frequent_kmers=[]
        
    for key, value in dict_kmer.items():
        if(value == max_value):
            list_frequent_kmers.append(key)
            
    return list_frequent_kmers

def main():
    with open('datasets/rosalind_ba1b.txt') as input_file:
        dna, k = input_file.read().strip().split('\n')
        
    list_frequent_kmers = frequent_words(dna, int(k))
    
    print(' '.join(list_frequent_kmers))
    
    with open('solutions/rosalind_ba1b.txt', 'w') as output_file:
        output_file.write(' '.join(list_frequent_kmers))
        
if(__name__=='__main__'):
    main()