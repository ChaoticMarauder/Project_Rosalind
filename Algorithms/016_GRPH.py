from rosalind import parse_fasta

def overlap_seq(seq1, seq2, k):
    
    if(seq1[len(seq1)-k:len(seq1)]==seq2[0:k]):
        return 1
    else:
        return 0
    
def overlap_graph(sequence_dict,k):
    
    overlap_list=[]
    
    for key1 in sequence_dict:
        for key2 in sequence_dict:
            if(sequence_dict[key1] != sequence_dict[key2]):
                edge=[]
                if(overlap_seq(sequence_dict[key1],sequence_dict[key2],k)==1):
                    
                    edge.append(key1)
                    edge.append(key2)
            
                    overlap_list.append(edge)
            
    return overlap_list

    
def main():
    
    file_name='datasets/rosalind_grph.txt'
    seq_dict = parse_fasta(file_name)
    
    k=3
    overlap_list=overlap_graph(seq_dict,k)
    
    for index in range(len(overlap_list)):
        print(' '.join(overlap_list[index]))
    
    with open('solutions/rosalind_grph.txt', 'w') as output_file:
        for index in range(len(overlap_list)):
            output_file.write(' '.join(overlap_list[index])+'\n')
        
if(__name__=='__main__'):
    main()