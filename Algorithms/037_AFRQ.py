import math

def hardy_weinberg(p_list):
    recessive_list=[]
    
    for idx in range(len(p_list)):
        r = math.sqrt(p_list[idx])
        h = 2-r
        
        recessive_allele_probability = r*h
        recessive_list.append(recessive_allele_probability)
        
    return recessive_list

def main():
    with open('datasets/rosalind_afrq.txt') as input_file:
        p_list = input_file.read().strip().split()
        
    p_list = list(map(float,p_list))
    
    recessive_list = hardy_weinberg(p_list)
    
    print(' '.join(list(map(str,recessive_list))))
    
    with open('solutions/rosalind_afrq.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str,recessive_list))))
         
if(__name__=="__main__"):
    main()