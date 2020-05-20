def summation(a,b):
    value=0
    for i in range(a,b+1):
        if(i%2!=0):
            value+=i
            
    return value

def main():
    
    with open('datasets/rosalind_ini4.txt') as input_file:
         a,b = input_file.read().strip().split()
     
    value=summation(int(a), int(b))
    print(str(value))
     
    with open('solutions/rosalind_ini4.txt', 'w') as output_file:
        output_file.write(str(value))
    
if(__name__=='__main__'):
    main()