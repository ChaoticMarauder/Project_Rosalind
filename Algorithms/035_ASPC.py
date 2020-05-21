from scipy.special import comb

def alternate_splice(n,m):
    alternate_splice_count=0
    for k in range(m, n+1):
        alternate_splice_count = (alternate_splice_count + comb(n, k, exact=True))%1000000
        
    return alternate_splice_count

def main():
    with open('datasets/rosalind_aspc.txt') as input_file:
        n, m = input_file.read().strip().split()
        
    alternate_splice_count=alternate_splice(int(n),int(m))
    
    print(str(alternate_splice_count))
    
    with open('solutions/rosalind_aspc.txt', 'w') as output_file:
        output_file.write(str(alternate_splice_count))
        
if(__name__=='__main__'):
    main()
