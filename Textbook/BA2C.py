from Chapter1 import NumberToPattern

def pattern_probability(pattern, profile):
    prob=1.0
    for i in range(len(pattern)):
        if(pattern[i]=='A'):
            prob=prob*profile[0][i]
        if(pattern[i]=='C'):
            prob=prob*profile[1][i]
        if(pattern[i]=='G'):
            prob=prob*profile[2][i]
        if(pattern[i]=='T'):
            prob=prob*profile[3][i]
            
    return prob

def profile_kmer(dna, k, profile):
    pattern_list=[]
    pattern_prob=[]
    
    for i in range(4**k):
        kmer = NumberToPattern(i,k)
        pattern_list.append(kmer)
        prob = pattern_probability(kmer, profile)
        pattern_prob.append(prob)
        
    pattern_prob, pattern_list = (list(t) for t in zip(*sorted(zip(pattern_prob, pattern_list))))

    for i in range(len(pattern_prob)-1,0,-1):
        for j in range(len(dna)-k+1):
            if(dna[j:j+k]==pattern_list[i]):
                return dna[j:j+k]
        
    
                
def main():
    with open('datasets/rosalind_ba2c.txt') as input_data:
         A = input_data.read().strip().split('\n')
        
    dna = A[0]
    
    k = int(A[1])
    
    profile = []
    
    profile = A[2:]
    
    for i in range(len(profile)):
        profile[i] = profile[i].split()
        profile[i] = list(map(float, profile[i]))
    
    
    motif = profile_kmer(dna, k, profile)
    
    print(motif)
    
    with open('solutions/rosalind_ba2c.txt', 'w') as output_file:
        output_file.write(motif)
        
if(__name__=='__main__'):
    main()