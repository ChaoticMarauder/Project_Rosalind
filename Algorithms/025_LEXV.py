from itertools import product

def lex_variable(list_letters, n):
    final_list_letters=['$'] + list_letters
    
    kmer_list=list(product(final_list_letters, repeat=n))
    kmer_lex_list=[]
    for kmer in kmer_list:
        
        if('$' not in kmer):
            kmer_lex_list.append(''.join(kmer))
        
        else:
            for i in range(1,n):
                if(''.join(kmer[i:n])=='$'*(n-i) and '$' not in kmer[:i]):
                    kmer_lex_list.append(''.join(kmer).replace('$',''))
    
    return kmer_lex_list

def main():
    with open('datasets/rosalind_lexv.txt') as input_file:
        list_letters, n = [line.strip() for line in input_file.readlines()]
        
    list_letters=list_letters.split()    
        
    kmer_lex_list=lex_variable(list_letters, int(n))    
    print('\n'.join(kmer_lex_list))
    
    
    with open('solutions/rosalind_lexv.txt', 'w') as output_file:
        output_file.write('\n'.join(kmer_lex_list))
        
if(__name__=='__main__'):
    main()