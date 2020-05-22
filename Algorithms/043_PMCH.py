from rosalind import parse_fasta
from rosalind import factorial

def perfect_matching(rna):
    a = rna.count('A')
    g = rna.count('G')
    
    num_perfect_matching = factorial(a)*factorial(g)
    
    return num_perfect_matching

def main():
    file_name='datasets/rosalind_pmch.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
        
    rna=seq_list[0]
    num_perfect_matching = perfect_matching(rna)
    
    print(str(num_perfect_matching))
    
    with open('solutions/rosalind_pmch.txt', 'w') as output_file:
        output_file.write(str(num_perfect_matching))
            
if(__name__=='__main__'):
    main()