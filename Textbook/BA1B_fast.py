from Chapter1 import PatternToNumber
from Chapter1 import NumberToPattern

def finding_frequent_patterns_sorting(dna, k):
    frequent_patterns=[]
    index=[]
    count=[]
    
    for i in range(len(dna)-k+1):
        pattern = dna[i:i+k]
        idx = PatternToNumber(pattern)
        cnt = 1
        index.append(idx)
        count.append(cnt)
        
    index.sort()
    
    for i in range(1,len(dna)-k+1):
        if(index[i]==index[i-1]):
            count[i]=count[i-1]+1
            
    max_count=max(count)
    for i in range(0,len(dna)-k+1):
        if(count[i]==max_count):
            pattern=NumberToPattern(index[i], k)
            frequent_patterns.append(pattern)
            
    return frequent_patterns

def main():
    with open('datasets/rosalind_ba1b.txt') as input_file:
        dna, k = input_file.read().strip().split('\n')
        
    list_frequent_kmers = finding_frequent_patterns_sorting(dna, int(k))
    
    print(' '.join(list_frequent_kmers))
    
    with open('solutions/rosalind_ba1b.txt', 'w') as output_file:
        output_file.write(' '.join(list_frequent_kmers))
        
if(__name__=='__main__'):
    main()
    
        