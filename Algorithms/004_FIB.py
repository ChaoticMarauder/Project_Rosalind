def rabbit_pairs(n,k):
    if(n==1):
        return 1
    if(n==2):
        return 1
    else:
        return k*rabbit_pairs(n-2,k)+rabbit_pairs(n-1,k)



def main():
    
    with open('datasets/rosalind_fib.txt') as input_data:
        n,k=map(int,input_data.read().strip().split())
    

#   n=5
#   k=3
    pairs=str(rabbit_pairs(n,k))
#    
#    
    print(pairs)
#    
    with open('solutions/rosalind_fib.txt', 'w') as output_data:
       output_data.write(pairs)
        
        
if (__name__=='__main__'):
    main()
