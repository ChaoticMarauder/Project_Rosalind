from itertools import product

def k_mer_lexicographic(bases, k):
    k_mers = [''.join(item) for item in product(bases, repeat=k)]
    return k_mers

def main():
    with open('datasets/rosalind_lexf.txt') as input_file:
        bases, k = input_file.readlines()
        
        bases=bases.split()
        k=int(k)
        
    k_mers=k_mer_lexicographic(bases, k)
    
    for k_mer in k_mers:
        print(k_mer+'\n')
            
    with open('solutions/rosalind_lexf.txt', 'w') as output_file:
        for k_mer in k_mers:
            output_file.write(k_mer+'\n')
            
            
if(__name__=='__main__'):
    main()