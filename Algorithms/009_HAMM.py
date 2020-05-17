def hamming_distance(s,t):
    count=0
    if(len(s)==len(t)):
        for i in range(len(s)):
            if(s[i]!=t[i]):
                count+=1
    else:
        return -1
        
    return count 

def main():
    with open('datasets/rosalind_hamm.txt') as input_data:
        s,t=input_data.read().strip().split('\n')
    
    ham_dist=hamming_distance(s,t)
    
    print(ham_dist)
    
    with open('solutions/rosalind_hamm.txt', 'w') as output_file:
        output_file.write(str(ham_dist))
        
if(__name__=='__main__'):
    main()    