def PatternCount(dna, pattern):
    count = 0
    len_pattern = len(pattern) 
    for i in range(len(dna)-len(pattern)):
        if(dna[i:i+len_pattern] == pattern):
            count+=1
    return count