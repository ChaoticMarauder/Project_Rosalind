from rosalind import parse_fasta

def check_motif(motif, seq_list):
    for seq in seq_list:
        if((motif not in seq) or (len(motif)>len(seq))):
            return False
    
    return True
        
def longest_shared_motif(seq_list):
    
    longest_motif=''
    length_first=len(seq_list[0])
    longest_motif_length=0
    
    for i in range(length_first):
        for j in range(length_first, i, -1):
            if(j-i<longest_motif_length):
                break
            elif(check_motif(seq_list[0][i:j],seq_list)):
                longest_motif=seq_list[0][i:j]
                longest_motif_length=j-i
    
    return longest_motif
    
    
def main():
    file_name='datasets/rosalind_lcsm.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
        
    shared_motif=longest_shared_motif(seq_list)
    print(shared_motif)
    
    with open('solutions/rosalind_lcsm.txt', 'w') as output_file:
        output_file.write(shared_motif)
    
    
if(__name__=='__main__'):
    main()    