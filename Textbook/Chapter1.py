def PatternCount(dna, pattern):
    count = 0
    len_pattern = len(pattern) 
    for i in range(len(dna)-len(pattern)):
        if(dna[i:i+len_pattern] == pattern):
            count+=1
    return count

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

def reverse_complement(dna):
    rc=[]
    for base in dna:
        if(base=='T'):
            rc.append('A')
        if(base=='A'):
            rc.append('T')
        if(base=='G'):
            rc.append('C')
        if(base=='C'):
            rc.append('G')    
    
    
    rc.reverse()
    return ''.join(rc)

def pattern_occurence(pattern, dna):
    list_occurences=[]
    k = len(pattern)
    for i in range(len(dna)-k):
        if(dna[i:i+k]==pattern):
            list_occurences.append(i)
            
    return list_occurences