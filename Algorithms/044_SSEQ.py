from rosalind import parse_fasta

def spliced_motif(dna, motif):
    splice_pos=[]
    i = 0
    
    for base in motif:
        
        while(base!=dna[i]):
            i = i+1
            
        i = i+1
        splice_pos.append(i)
        
    return splice_pos

def main():        
    file_name='datasets/rosalind_sseq.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
    
    dna = seq_list[0]
    motif = seq_list[1]
    
    splice_pos = spliced_motif(dna, motif)
    
    print(' '.join(map(str,(splice_pos))))
    
    with open('solutions/rosalind_sseq.txt', 'w') as output_file:
        output_file.write(' '.join(map(str,(splice_pos))))
        
if(__name__=='__main__'):
    main()