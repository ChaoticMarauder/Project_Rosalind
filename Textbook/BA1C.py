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

def main():
    with open('datasets/rosalind_ba1c.txt') as input_file:
        dna = input_file.read().strip()
        
    rc = reverse_complement(dna)
    
    print(rc)
    
    with open('solutions/rosalind_ba1c.txt', 'w') as output_file:
        output_file.write(rc)
        
if(__name__=='__main__'):
    main()