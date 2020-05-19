from rosalind import translate
from rosalind import transcribe
from rosalind  import reverse_complement
from rosalind import parse_fasta

def last_stop_codon(rna):
    
    last_stop=0
    for i in range(0,len(rna)-3,3):
        codon=rna[i:i+3]
        if(codon=="UAA" or codon=="UAG" or codon=="UGA"):
            last_stop=i+2
            continue
            
    return last_stop
        
def orf_translate(dna):
    
    rna=''.join(transcribe(dna))
    reverse_rna=''.join(transcribe(''.join(reverse_complement(dna))))
    
    orf_list=[]
    
    last_stop_rna=[]
    last_stop_rna.append(last_stop_codon(rna))
    last_stop_rna.append(last_stop_codon(rna[1:]))
    last_stop_rna.append(last_stop_codon(rna[2:]))
    
    last_stop_reverse=[]
    last_stop_reverse.append(last_stop_codon(reverse_rna))
    last_stop_reverse.append(last_stop_codon(reverse_rna[1:]))
    last_stop_reverse.append(last_stop_codon(reverse_rna[2:]))
    
    first_strand=[]
    second_strand=[]
    
    for i in range(3):
        first_strand.append(rna[i:last_stop_rna[i]+i+1])
        second_strand.append(reverse_rna[i:last_stop_reverse[i]+i+1])
    
    print(first_strand)
    print(second_strand)
    
    for j in range(3):
        for k in range(0,len(first_strand[j]),3):
    
            seq1=translate(first_strand[j][k:len(first_strand[j])-3])
            if(len(seq1)!=0):
                orf_list.append(seq1)
                print(seq1)
                print(j)
                print(k)
        
        for k in range(0,len(second_strand[j]),3):    
        
            seq2=translate(second_strand[j][k:len(second_strand[j])-3])
            if(len(seq2)!=0):
                orf_list.append(seq2)
            
    
    
    
    possible_orfs = [] 
    
    for seq in orf_list: 
        if seq not in possible_orfs:
            possible_orfs.append(seq) 
            
            
    return possible_orfs


def main():
    file_name='datasets/rosalind_orf.txt'
    seq_dict = parse_fasta(file_name)
    
    for key in seq_dict:
        dna=seq_dict[key]
    

    orfs=orf_translate(dna)
    
    print('\n'.join(orfs))
    
    with open('solutions/rosalind_orf.txt', 'w') as output_file:
        output_file.write('\n'.join(orfs))
        
if(__name__=='__main__'):
    main()
    