def mortalfibonacciRabbits(n, m):
    rabbits = [0, 1, 1]
    for i in range(3, n + 1):
        if(i <= m):
            total = rabbits[i - 1] + rabbits[i - 2]
        elif(i == m + 1):
            total = rabbits[i - 1] + rabbits[i - 2] - 1
        else:
            total = rabbits[i - 1] + rabbits[i - 2] - rabbits[i - m - 1]
        rabbits.append(total)
    return (rabbits[n])
   



def main():
    
    with open('datasets/rosalind_fibd.txt') as input_data:
        n,m=map(int,input_data.read().strip().split())
    

    pairs=str(mortalfibonacciRabbits(n, m))
   
    print(pairs)
    
    with open('solutions/rosalind_fibd.txt', 'w') as output_data:
       output_data.write(pairs)
        
        
if (__name__=='__main__'):
    main()
