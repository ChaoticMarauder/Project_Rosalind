import numpy as np
from rosalind import parse_fasta

def consensus_profile(seq_list):
    seq_length=len(seq_list[0])
    consensus_array=np.zeros([4,seq_length])
    
    for i in range(len(seq_list)):
        for j in range(seq_length):
            if(seq_list[i][j]=='A'):
                consensus_array[0][j]+=1
            if(seq_list[i][j]=='C'):
                consensus_array[1][j]+=1
            if(seq_list[i][j]=='G'):
                consensus_array[2][j]+=1
            if(seq_list[i][j]=='T'):
                consensus_array[3][j]+=1
                
    consensus_array=consensus_array.astype(int)
                
    dict_consensus={'A':consensus_array[0], 'C':consensus_array[1],
                        'G':consensus_array[2], 'T':consensus_array[3] }
        
    return dict_consensus, consensus_array


def base_mapping(n):
    if(n==0):
        base='A'
    if(n==1):
        base='C'
    if(n==2):
        base='G'
    if(n==3):
        base='T'
        
    return base
    
def consensus_sequence(array_cons):
    
    consensus_seq=""
    
    for i in range(np.size(array_cons,1)):
        max_value=0
        index=0
        for j in range(np.size(array_cons,0)):
            if(array_cons[j][i]>max_value):
                max_value=array_cons[j][i]
                index=j
            
        consensus_seq=consensus_seq+base_mapping(index)
            
    return consensus_seq
        
    
def main():
    file_name='datasets/rosalind_cons.txt'
    
    seq_dict = parse_fasta(file_name)
    
    sequence_list=[]
    for key in seq_dict:
        sequence_list.append(seq_dict[key])
    
    dict_cons, array_cons=consensus_profile(sequence_list)
    
    con_seq=consensus_sequence(array_cons)
    
    print(con_seq)
    for key in dict_cons:
        print(key+': '+' '.join(map(str,(dict_cons[key]))))
   
    
    
    with open('solutions/rosalind_cons.txt', 'w') as output_file:
        output_file.write(con_seq+'\n')
        for key in dict_cons:
            output_file.write(key+': '+' '.join(map(str,(dict_cons[key])))+'\n')
        
    
if(__name__=="__main__"):
    main()