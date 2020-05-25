from Chapter1 import reverse_complement

def transcribe(dna):
    rna=[]
    for base in dna:
        if(base=='T'):
            rna.append('U')
        else:
            rna.append(base)
    
    return ''.join(rna)

def encode_peptide(dna, peptide):
    
    dict_codon = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"X", "UAG":"X",
       "UGU":"C", "UGC":"C", "UGA":"X", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
    
    rna = transcribe(dna)
    rev_comp = reverse_complement(dna)
    rev_comp_rna = transcribe(rev_comp)
    
    prot_candidates=['','','','','','']
    for i in range(int(len(rna)/3)):
        if(len(rna[3*i:3*i+3])==3):
            prot_candidates[0]=prot_candidates[0]+dict_codon[rna[3*i:3*i+3]]
        if(len(rna[3*i+1:3*i+4])==3):
            prot_candidates[1]=prot_candidates[1]+dict_codon[rna[3*i+1:3*i+4]]
        if(len(rna[3*i+2:3*i+5])==3):
            prot_candidates[2]=prot_candidates[2]+dict_codon[rna[3*i+2:3*i+5]]
        if(len(rev_comp_rna[3*i:3*i+3])==3):
            prot_candidates[3]=prot_candidates[3]+dict_codon[rev_comp_rna[3*i:3*i+3]]
        if(len(rev_comp_rna[3*i+1:3*i+4])==3):
            prot_candidates[4]=prot_candidates[4]+dict_codon[rev_comp_rna[3*i+1:3*i+4]]
        if(len(rev_comp_rna[3*i+2:3*i+5])==3):
            prot_candidates[5]=prot_candidates[5]+dict_codon[rev_comp_rna[3*i+2:3*i+5]]
        
    list_dna=[]    
    k = len(peptide)    
    for i in range(0,3):
        for j in range(len(prot_candidates[i])-k+1):
            if(prot_candidates[i][j:j+k] == peptide):
                pattern = dna[3*j+i:3*(j+k-1)+i+3]
                list_dna.append(pattern)
    
    for i in range(3,6):
        for j in range(len(prot_candidates[i])-k+1):
            if(prot_candidates[i][j:j+k] == peptide):
                pattern = rev_comp[3*j+i-3:3*(j+k-1)+i-3+3]
                list_dna.append(reverse_complement(pattern))       

    return list_dna


def main():
    with open('datasets/rosalind_ba4b.txt') as input_file:
        dna, peptide = input_file.read().strip().split('\n')
        
    list_dna = encode_peptide(dna, peptide)
    
    print('\n'.join(list_dna))
    
    with open('solutions/rosalind_ba4b.txt', 'w') as output_file:
        output_file.write('\n'.join(list_dna))
        
if(__name__=='__main__'):
    main()
        
    
    