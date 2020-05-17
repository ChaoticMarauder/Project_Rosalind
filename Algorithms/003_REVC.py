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
    return rc




def main():
    
    with open('datasets/rosalind_revc.txt') as input_data:
        dna=input_data.read().strip().split()
    
#   dna='AAAACCCGGT'
    rc=reverse_complement(dna)
    
    
    print(''.join(rc))
    
    with open('solutions/rosalind_revc.txt', 'w') as output_data:
        output_data.write(''.join(map(str,(rc))))
        
        
if (__name__=='__main__'):
    main()
