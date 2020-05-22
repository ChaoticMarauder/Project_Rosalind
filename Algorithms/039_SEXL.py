def sex_linked_inheritance(allele_frequency_list):
    carrier_list=[]
    for idx in range(len(allele_frequency_list)):
        carrier = round(2*allele_frequency_list[idx]*(1-allele_frequency_list[idx]),3)
        carrier_list.append(carrier)
    
    return carrier_list

def main():
    with open('datasets/rosalind_sexl.txt') as input_file:
        allele_frequency_list = input_file.read().strip().split()
        
    allele_frequency_list=list(map(float,allele_frequency_list))
      
    carrier_list = sex_linked_inheritance(allele_frequency_list)  
    print(' '.join(list(map(str,(carrier_list)))))
    
    with open('solutions/rosalind_sexl.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str,(carrier_list)))))
        
if(__name__=='__main__'):
    main()