from rosalind import parse_fasta
from rosalind import reverse_complement

def reverse_palindrome(dna):
    if(dna==''.join(reverse_complement(dna))):
        return 1
    else:
        return 0
    
def restriction_sites(dna):
    len_dna=len(dna)
    
    site_list=[]
    length_list=[]
    
    for j in range(4,13):
        for i in range(len_dna-j+1):
                if(reverse_palindrome(dna[i:i+j])==1):
                    site_list.append(i+1)
                    length_list.append(j)
                    
    return site_list, length_list

def main():
    
    file_name='datasets/rosalind_revp.txt'
    seq_dict = parse_fasta(file_name)
    
    for key in seq_dict:
        dna=seq_dict[key]
        
    restrict_sites, len_sites=restriction_sites(dna)
    
    for i in range(len(restrict_sites)):
         print(str(restrict_sites[i])+' '+str(len_sites[i]))
         
    with open('solutions/rosalind_revp.txt', 'w') as output_file:
        for i in range(len(restrict_sites)):
            output_file.write(str(restrict_sites[i])+' '+str(len_sites[i])+'\n')
    
if(__name__=='__main__'):
    main()