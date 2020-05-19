from itertools import permutations

def factorial(n):
    if(n<=1):
        return 1
    else:
        return n*factorial(n-1)
    
def rearrangements(n):
    num_set=range(1,n+1)
    permutation_list = list(permutations(num_set))
    
    return permutation_list

def main():
    
    with open('datasets/rosalind_perm.txt') as input_file:
        n=int(input_file.read())
        
    list_perm=rearrangements(n)
    
    print(factorial(n))
    for idx in range(len(list_perm)):
        print(' '.join(map(str, (list_perm[idx]))))
    
    with open('solutions/rosalind_perm.txt', 'w') as output_file:
        output_file.write(str(factorial(n))+'\n')
        for idx in range(len(list_perm)):
             output_file.write(' '.join(map(str, (list_perm[idx])))+'\n')
if(__name__=='__main__'):
    main()
