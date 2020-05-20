from scipy.special import comb

def independent_alleles(k, N):
    prob_N=0.0
    total_offspring=2**k
    
    for i in range(N, total_offspring):
        #calculating the success of at least N in among 2**k trials-->binomial distribution
        prob_N=prob_N+comb(total_offspring,i)*((0.25)**i)*((0.75)**(total_offspring-i))
    
    prob_N=round(prob_N,4)
    return prob_N

def main():
    with open('datasets/rosalind_lia.txt') as input_file:
        k, N = input_file.read().strip().split()
        
    
    prob_N=independent_alleles(int(k),int(N))
    
    print(prob_N)
    
    with open('solutions/rosalind_lia.txt', 'w') as output_file:
        output_file.write(str(prob_N))
        
if(__name__=='__main__'):
    main()