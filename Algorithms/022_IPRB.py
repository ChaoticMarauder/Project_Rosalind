from scipy.special import comb

def mendel_first_law(dominant, heterozygous, recessive):
    num_offspring = comb(dominant+heterozygous+recessive, 2)
    
    num_recessive_offspring = (comb(recessive, 2) + 0.5*recessive*heterozygous
                              + 0.25*comb(heterozygous, 2))
    
    num_dominant_offspring = (comb(dominant, 2) + 0.75*comb(heterozygous, 2) 
                             + dominant*heterozygous + dominant*recessive 
                             + 0.5*recessive*heterozygous)
    
    ratio_recessive = num_recessive_offspring/num_offspring
    
    ratio_dominant = 1-ratio_recessive
    ratio_dominant = num_dominant_offspring/num_offspring
    
    return ratio_dominant

def main():
    with open('datasets/rosalind_iprb.txt') as input_file:
        k,m,n = input_file.read().strip().split()
        
    ratio_dominant=mendel_first_law(int(k),int(m),int(n))
    
    print(ratio_dominant)
    
    with open('solutions/rosalind_iprb.txt', 'w') as output_file:
        output_file.write(str(ratio_dominant))
        
if(__name__=='__main__'):
    main()
        
        