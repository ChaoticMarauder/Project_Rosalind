def dna_base_count(dna):
    count_list=[]
    count_A=0
    count_C=0
    count_G=0
    count_T=0
    for base in dna:
        if(base=='A'):
            count_A+=1
        if(base=='C'):
            count_C+=1
        if(base=='G'):
            count_G+=1
        if(base=='T'):
            count_T+=1
    
    count_list=[count_A, count_C, count_G, count_T]
    
    return count_list
 

def main():
    
    with open('datasets/rosalind_dna.txt') as input_data:
        dna=input_data.read().strip()
    
    count=dna_base_count(dna)
    
    print(' '.join(map(str,(count))))
    
    with open('solutions/rosalind_dna.txt', 'w') as output_data:
        output_data.write(' '.join(map(str,(count))))
        
        
if (__name__=='__main__'):
    main()
