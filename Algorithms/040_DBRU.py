from rosalind import reverse_complement

def de_bruijn_graph(seq_list):
    
    final_seq_list=[]
    
    for seq in seq_list:
        if seq not in final_seq_list:
            final_seq_list.append(seq)
    
    reverse_complement_seq_list=[]
    
    for seq in final_seq_list:
        reverse_complement_seq_list.append(''.join(reverse_complement(seq)))
        
    adjacency_list=[]
    k=len(final_seq_list[0])
    
    for seq in final_seq_list:
        edge=['','']
        edge[0]=seq[0:k-1]
        edge[1]=seq[1:k]
        if edge not in adjacency_list:
            adjacency_list.append(edge)
    
    for seq in reverse_complement_seq_list:
        edge=['','']
        edge[0]=seq[0:k-1]
        edge[1]=seq[1:k]
        if edge not in adjacency_list:
            adjacency_list.append(edge)
            
    return adjacency_list

def main():
    with open('datasets/rosalind_dbru.txt') as input_file:
        seq_list = input_file.read().strip().split('\n')
        
    adjacency_list=de_bruijn_graph(seq_list)
      
    for edge in adjacency_list:
        print('('+', '.join(edge)+')')
    
    with open('solutions/rosalind_dbru.txt', 'w') as output_file:
        for edge in adjacency_list:
            output_file.write('('+', '.join(edge)+')'+'\n')
            
if(__name__=='__main__'):
    main()