def NumberToSymbol(index):
    base=''
    if(index==0):
        base = 'A'
    if(index==1):
        base = 'C'
    if(index==2):
        base = 'G'
    if(index==3):
        base = 'T'
        
    return base

def NumberToPattern(index, k):
       
    if(k==1):
        return NumberToSymbol(index)
    
    r=index % 4
    prefix_index=int(index / 4)
    
    return NumberToPattern(prefix_index, k-1) + NumberToSymbol(r) 
    
    
def main():
    with open('datasets/rosalind_ba1m.txt') as input_file:
        index, k = input_file.read().strip().split('\n')
        
    index=int(index)
    k=int(k)
    
    print(index)
    print(k)
    pattern = NumberToPattern(index, k)
    
    print(pattern)
    
    with open('solutions/rosalind_ba1m.txt', 'w') as output_file:
        output_file.write(pattern)

if(__name__=='__main__'):
    main()

