def matching_random_motif(N, gc_content, seq):
    GC_count=0
    AT_count=0
    
    for i in range(len(seq)):
        if(seq[i]=='G' or seq[i]=='C'):
            GC_count+=1
        else:
            AT_count+=1
            
    prob_seq = (gc_content/2)**GC_count*((1-gc_content)/2)**AT_count
    prob_not_motif = (1-prob_seq)**N
    
    prob_motif = 1-prob_not_motif
    
    return prob_motif

def main():
    with open('datasets/rosalind_rstr.txt') as input_file:
        num, seq = [line.strip() for line in input_file.readlines()]
        
    num = num.split()
    
    N = int(num[0])
    gc_content = float(num[1])
    
    prob = matching_random_motif(N, gc_content, seq)
    
    print(prob)
    
    with open('solutions/rosalind_rstr.txt', 'w') as output_file:
        output_file.write(str(prob)+' ')
        
if(__name__=='__main__'):
    main() 