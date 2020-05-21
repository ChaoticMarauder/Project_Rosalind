def mrna_to_protein(protein):
    codon_dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
       "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
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
    
    mrna_possibilities = list(codon_dict.values()).count("STOP")

    for aa in protein:
        mrna_possibilities = (mrna_possibilities*list(codon_dict.values()).count(aa))% 1000000

    return mrna_possibilities 

def main():
    with open('datasets/rosalind_mrna.txt') as input_file:
        protein=input_file.read().strip()
        
    mrna_possibilities=mrna_to_protein(protein)
    
    print(str(mrna_possibilities))
    
    with open('solutions/rosalind_mrna.txt', 'w') as output_file:
        output_file.write(str(mrna_possibilities))
        
if(__name__=='__main__'):
    main() 