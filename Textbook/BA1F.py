def gc_genome_skew(dna):
    gc_skew=[]
    gc_diff=0
    
    for i in range(len(dna)):
        gc_skew.append(gc_diff)
        if dna[i]=='C':
            gc_diff-=1
        if dna[i]=='G':
            gc_diff+=1
    
    oric_list=[]
    min_value = min(gc_skew)
    
    for i in range(len(gc_skew)):
        if(gc_skew[i]==min_value):
            oric_list.append(i)
    
    return oric_list

def main():
    with open('datasets/rosalind_ba1f.txt') as input_file:
        dna = input_file.read().strip()
        
    oric_list = gc_genome_skew(dna)
    
    print(' '.join(list(map(str, oric_list))))
    
    with open('solutions/rosalind_ba1f.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str, oric_list))))
        
if(__name__=='__main__'):
    main()