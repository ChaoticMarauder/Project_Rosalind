from rosalind import parse_fasta

def max_overlap_sequence(seq_list):
    
    max_overlap=-1
    length_list=len(seq_list)
    
    for first_id in range(length_list):
        for second_id in range(length_list):
            if(first_id==second_id):
                continue
            else:
                first_string, second_string = seq_list[first_id], seq_list[second_id] 
                
                k=0
                
                while(first_string[k:]!=second_string[0:len(first_string[k:])]):
                    k=k+1
                
                if(len(first_string)-k > max_overlap):
                    max_overlap=len(first_string)-k
                    max_index=[first_id, second_id]
                    
    return [seq_list[max_index[0]]+seq_list[max_index[1]][max_overlap:]]+[seq_list[i] for i in range(length_list) if i not in max_index]

def shortest_super_sequence(seq_list):
    
    while(len(seq_list)>1):
        seq_list=max_overlap_sequence(seq_list)
        
    return seq_list[0]

def main():
    
    file_name='datasets/rosalind_long.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
        
    final_sequence=shortest_super_sequence(seq_list)
    
    print(final_sequence)
    print(len(final_sequence))
    
    with open('solutions/rosalind_long.txt', 'w') as output_file:
        output_file.write(final_sequence)
        
if(__name__=='__main__'):
    main()
    