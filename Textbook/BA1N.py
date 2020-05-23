from Chapter1 import hamming_distance

def neighbours(pattern, d):
    if d==0:
        return [pattern]
    if len(pattern)==1:
        return ['A','C','G','T']
    
    neighbourhood=[]
    
    suffix_neighbours = neighbours(pattern[1:len(pattern)], d)
    
    for pattern_suffix in suffix_neighbours:
        if hamming_distance(pattern[1:len(pattern)],pattern_suffix) < d:
            for nucleotide in ['A','C','G','T']:
                pattern_new = nucleotide + pattern_suffix
                neighbourhood.append(pattern_new)
        else:
            pattern_new = pattern[0] + pattern_suffix
            neighbourhood.append(pattern_new)
            
    return neighbourhood

def main():
    with open('datasets/rosalind_ba1n.txt') as input_data:
        pattern, d = input_data.read().strip().split('\n')
        
    d = int(d)
    neighbourhood = neighbours(pattern, d)
    
    print('\n'.join(neighbourhood))
    
    with open('solutions/rosalind_ba1n.txt', 'w') as output_file:
        output_file.write('\n'.join(neighbourhood))
        
if(__name__=='__main__'):
    main()