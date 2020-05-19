from itertools import product, permutations

def signed_permutations(n):
    signed_ints_list=[]
    
    for i in range(1,n+1):
        signed_ints_list.append([i,-1*i])
    
    product_list=map(list,list(product(*signed_ints_list)))
    
    signed_permutations_list=[]
    for pro in product_list:
        signed_permutations_list+=map(list,list(permutations(pro)))
        
    return signed_permutations_list

def main():
    
    with open('datasets/rosalind_sign.txt') as input_file:
        n=input_file.read()
        
    permutation_list=signed_permutations(int(n))
    
    print(len(permutation_list))
    for item in permutation_list:
        print(' '.join(map(str,item)))
    
    with open('solutions/rosalind_sign.txt', 'w') as output_file:
        output_file.write(str(len(permutation_list))+'\n')
        for item in permutation_list:
            output_file.write(' '.join(map(str,item))+'\n')
        
if(__name__=='__main__'):
    main()