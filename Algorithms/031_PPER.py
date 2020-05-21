def partial_perm(n, k):
    num_perm=1
    for i in range(n, n-k, -1):
        num_perm = (num_perm*i)%1000000
        
    return num_perm

def main():
    with open('datasets/rosalind_pper.txt') as input_file:
        n, k = input_file.read().strip().split()
        
    num_perm = partial_perm(int(n), int(k))
    
    print(str(num_perm))
    
    with open('solutions/rosalind_pper.txt', 'w') as output_file:
       output_file.write(str(num_perm))
       
if(__name__=='__main__'):
    main()