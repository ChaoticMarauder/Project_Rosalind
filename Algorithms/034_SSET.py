def num_subsets(n):
    sset=1
    for i in range(n):
        sset=(sset*2)%1000000
    
    return sset

def main():
    with open('datasets/rosalind_sset.txt') as input_file:
        n = int(input_file.read().strip())
        
    sset=num_subsets(n)
    
    print(str(sset))
    
    with open('solutions/rosalind_sset.txt', 'w') as output_file:
        output_file.write(str(sset))
        
if(__name__=='__main__'):
    main()