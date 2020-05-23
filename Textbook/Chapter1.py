
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

#faster and more efficient version of the frequent words function
def finding_frequent_patterns_sorting(dna, k):
    frequent_patterns=[]
    index=[]
    count=[]
    
    for i in range(len(dna)-k):
        pattern = dna[i:i+k]
        idx = PatternToNumber(pattern)
        cnt = 1
        index.append(idx)
        count.append(cnt)
        
    index.sort()
    
    for i in range(1,len(dna)-k):
        if(index[i]==index[i-1]):
            count[i]=count[i-1]+1
            
    max_count=max(count)
    for i in range(0,len(dna)-k):
        if(count[i]==max_count):
            pattern=NumberToPattern(index[i], k)
            frequent_patterns.append(pattern)
            
    return frequent_patterns

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

def SymbolToNumber(base):
    if(base=='A'):
        return 0
    if(base=='C'):
        return 1
    if(base=='G'):
        return 2
    if(base=='T'):
        return 3
    
def PatternToNumber(pattern):
    if(len(pattern)==0):
        return 0
    
    last_symbol = pattern[-1]
    prefix = pattern[0:len(pattern)-1]
    
    return 4*PatternToNumber(prefix)+SymbolToNumber(last_symbol)

def NumberToSymbol(index):
    base=''
    if(index==0):
        base = 'A'
    if(index==1):
        base = 'C'
    if(index==2):
        base = 'G'
    if(index==3):
        base = 'T'
        
    return base

def NumberToPattern(index, k):
       
    if(k==1):
        return NumberToSymbol(index)
    
    r=index % 4
    prefix_index=int(index / 4)
    
    return NumberToPattern(prefix_index, k-1) + NumberToSymbol(r) 


def computing_frequencies(dna, k):
    frequency_array=[]
    for i in range(4**k-1):
        frequency_array.append(0)
        
    for i in range(len(dna)-k):
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