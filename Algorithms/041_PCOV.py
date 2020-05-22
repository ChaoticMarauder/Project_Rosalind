def cyclic_assembly(seq_list):
    final_seq_list=[]
    
    for seq in seq_list:
        if seq not in final_seq_list:
            final_seq_list.append(seq)
    
    adjacency_list=[]
    k=len(final_seq_list[0])
    
    for seq in final_seq_list:
        edge=['','']
        edge[0]=seq[0:k-1]
        edge[1]=seq[1:k]
        if edge not in adjacency_list:
            adjacency_list.append(edge)
            
    temp_kmer = adjacency_list.pop(0)
    cyclic_seq = temp_kmer[0][-1]
    while(len(adjacency_list)!=0):
        cyclic_seq = cyclic_seq + temp_kmer[1][-1]
        [idx]=[index for index, edges in enumerate(adjacency_list) if(edges[0]==temp_kmer[1])]
        temp_kmer = adjacency_list.pop(idx)
        
    return cyclic_seq

def main():
    with open('datasets/rosalind_pcov.txt') as input_file:
        seq_list = input_file.read().strip().split('\n')
        
    cyclic_seq=cyclic_assembly(seq_list)
      
    print(cyclic_seq)
    
    with open('solutions/rosalind_pcov.txt', 'w') as output_file:
        output_file.write(cyclic_seq)
            
if(__name__=='__main__'):
    main()