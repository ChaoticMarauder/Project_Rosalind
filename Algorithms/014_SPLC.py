from rosalind import parse_fasta
from rosalind import translate
from rosalind import transcribe

def rna_splicing(dna, intron_list):
    
    for intron in intron_list:
        dna=dna.replace(intron,'')
    
    spliced_rna = ''.join(transcribe(dna))
    
    return spliced_rna

def main():
    
    file_name='datasets/rosalind_splc.txt'
    seq_dict = parse_fasta(file_name)
    
    list_seq=[]
    
    for key in seq_dict:
        list_seq.append(seq_dict[key])
        
    dna=list_seq[0]
    intron_list=list_seq[1:]
    
    spliced_rna=rna_splicing(dna, intron_list)
    
    protein, start_codon, end_codon=translate(spliced_rna)
    
    print(protein)
    
    with open('solutions/rosalind_splc.txt', 'w') as output_file:
        output_file.write(protein)
        
if(__name__=='__main__'):
    main()