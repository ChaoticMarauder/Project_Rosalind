def counting_peptides(mass):
    dict_weights={'A':71,'C':103,'D':115,'E':129,'F':147,'G':57,'H':137,
                  'I':113,'K':128,'L':113,'M':131,'N':114,'P':97,'Q':128,
                  'R':156,'S':87,'T':101,'V':99,'W':186,'Y':163}
    
    list_weights=[]
    
    for key, value in dict_weights.items():
        if value not in list_weights:
            list_weights.append(value)
            
            
    len_list = len(list_weights) 
    
    count_peptides=[]
    
    for i in range(mass+1):
        count_peptides.append(0)
    
    count_peptides[0] = 1

    for i in range(1, mass+1):
        for j in range(len_list):
            if(i >= list_weights[j]):
                count_peptides[i] = count_peptides[i] + count_peptides[i-list_weights[j]]
                
    return count_peptides[mass]

def main():
    with open('datasets/rosalind_ba4d.txt') as input_file:
        mass = int(input_file.read().strip())
        
    count_peptides = counting_peptides(mass)
    
    print(str(count_peptides))
    
    with open('solutions/rosalind_ba4d.txt', 'w') as output_file:
        output_file.write(str(count_peptides))
        
if(__name__=='__main__'):
    main()