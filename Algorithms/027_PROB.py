import math

def prob_gc_seq(seq, gc_content):
    prob=0.0
    for i in range(len(seq)):
        
        if(seq[i]=='G' or seq[i]=='C'):
            prob=prob+math.log10((gc_content/2))
        
        if(seq[i]=='A' or seq[i]=='T'):
            prob=prob+math.log10((1-gc_content)/2)
            
    log_prob=round(prob,4)
    
    return log_prob

def list_prob_gc_seq(seq, list_gc_content):
    list_log_prob=[]
    for i in range(len(list_gc_content)):
        gc_content=list_gc_content[i]
        list_log_prob.append(prob_gc_seq(seq,gc_content))
        
    return list_log_prob

def main():
    
    with open('datasets/rosalind_prob.txt') as input_file:
        seq, list_gc_content = [line.strip() for line in input_file.readlines()]
    
    list_gc_content = list_gc_content.split()
    list_gc_content = list(map(float,(list_gc_content)))
    
    list_log_prob =list_prob_gc_seq(seq, list_gc_content)
    
    for prob in list_log_prob:
        print(str(prob))
    
    with open('solutions/rosalind_prob.txt', 'w') as output_file:
        for prob in list_log_prob:
            output_file.write(str(prob)+' ')
        
if(__name__=='__main__'):
    main() 
    