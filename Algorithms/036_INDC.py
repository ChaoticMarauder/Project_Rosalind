from scipy.special import comb
import math

def independent_segregation(n):
    
    probability = 2**(2*n)
    
    log_prob_array = [-2*n*math.log10(2)]*2*n
    
    for k in range(2*n):
        probability = probability - comb(2*n, k, exact=True)
        log_prob_array[k]=log_prob_array[k] + math.log10(probability)
        
    return log_prob_array

def main():
    with open('datasets/rosalind_indc.txt') as input_file:
        n = input_file.read().strip()
        
    n = int(n)
      
    log_prob_list = independent_segregation(n)  
    print(' '.join(list(map(str,(log_prob_list)))))
    
    with open('solutions/rosalind_indc.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str,(log_prob_list)))))
        
if(__name__=='__main__'):
    main()