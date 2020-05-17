def transcribe(dna):
    rna=[]
    for base in dna:
        if(base=='T'):
            rna.append('U')
        else:
            rna.append(base)
    
    return rna




def main():
    
    with open('datasets/rosalind_rna.txt') as input_data:
        dna=input_data.read().strip()
    
#    dna='GATGGAACTTGACTACGTAAATT'
    rna=transcribe(dna)
    
    
    print(''.join(rna))
    
    with open('solutions/rosalind_rna.txt', 'w') as output_data:
        output_data.write(''.join(map(str,(rna))))
        
        
if (__name__=='__main__'):
    main()
