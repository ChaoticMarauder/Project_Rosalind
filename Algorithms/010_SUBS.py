def dna_motif(dna, motif):
    m_len=len(motif)
    list_motifs=[]
    for i in range(len(dna)-m_len):
        if(dna[i:i+m_len]==motif):
            list_motifs.append(i+1)
            
    return list_motifs

def main():
    
    with open('datasets/rosalind_subs.txt') as input_data:
        dna,motif=input_data.read().strip().split('\n')
    
#    dna='GATATATGCATATACTT'
#    motif='ATAT'
    motifs=dna_motif(dna, motif)
    
    
    print(' '.join(map(str,(motifs))))
    
    with open('solutions/rosalind_subs.txt', 'w') as output_data:
        output_data.write(' '.join(map(str,(motifs))))
    
if(__name__=='__main__'):
    main()    