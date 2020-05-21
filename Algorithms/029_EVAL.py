def motif_prob(seq, gc_content):
    GC_count=0
    AT_count=0
    
    for i in range(len(seq)):
        if(seq[i]=='G' or seq[i]=='C'):
            GC_count+=1
        else:
            AT_count+=1
            
    prob_seq = (gc_content/2)**GC_count*((1-gc_content)/2)**AT_count
    
    return prob_seq
    
def expected_restriction_sites(N, seq, gc_content):
    #seq is the restriction site sequence
    n = N - len(seq) + 1
    motif_probability = motif_prob(seq, gc_content)
    expected_sites = n*motif_probability
    
    return expected_sites

def main():
    with open('datasets/rosalind_eval.txt') as input_file:
        N, seq, gc_content_list  = [line.strip() for line in input_file.readlines()]
        
    N = int(N)
    gc_content_list = gc_content_list.split()
    gc_content_list = list(map(float,gc_content_list))
    
    expected_sites_list=[]    
    for i in range(len(gc_content_list)):
        expected_sites = expected_restriction_sites(N, seq, gc_content_list[i])
        expected_sites_list.append(expected_sites)
        
    print(' '.join(list(map(str,(expected_sites_list)))))
    
    with open('solutions/rosalind_eval.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str,(expected_sites_list)))))
        
if(__name__=='__main__'):
    main() 