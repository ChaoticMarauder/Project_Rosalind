from rosalind import k_mer_lexicographic
from rosalind import parse_fasta
import numpy as np

def k_mer_sequence(dna):
    kmer_list=k_mer_lexicographic('ACGT',4)
    kmer_count=np.zeros(len(kmer_list))
    
    for i in range(len(dna)-3):
        kmer=dna[i:i+4]
        kmer_count[kmer_list.index(kmer)]+=1
    
    kmer_count=kmer_count.astype(int)
    return kmer_count

def main():
    file_name='datasets/rosalind_kmer.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
    
    dna=seq_list[0]
    
    kmer_count=k_mer_sequence(dna)
    
    print(' '.join(map(str,(kmer_count))))
    
    with open('solutions/rosalind_kmer.txt', 'w') as output_file:
        output_file.write(' '.join(map(str,(kmer_count))))
        
if(__name__=='__main__'):
    main()