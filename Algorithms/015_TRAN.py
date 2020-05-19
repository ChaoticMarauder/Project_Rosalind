from rosalind import parse_fasta

def transit_transverse_ratio(dna_list):
    
    transition=0
    transversion=0
    
    length=len(dna_list[0])
    
    for i in range(length):
        A=dna_list[0][i]
        B=dna_list[1][i]
        
        if(A!=B):
            if((A=='A' and B=='G') or (A=='G' and B=='A')):
                transition+=1
            elif((A=='C' and B=='T') or (A=='T' and B=='C')):
                transition+=1
            else:
                transversion+=1 
                
    if(transition!=0 and transversion!=0):
        ratio=round(transition/transversion,11)
        return ratio
    else:
        return 0

def main():
    file_name='datasets/rosalind_tran.txt'
    seq_dict = parse_fasta(file_name)
    
    list_seq=[]
    
    for key in seq_dict:
        list_seq.append(seq_dict[key])
        
        
    ratio=transit_transverse_ratio(list_seq)
    
    print(ratio)
    
    with open('solutions/rosalind_tran.txt', 'w') as output_file:
        output_file.write(str(ratio))
    
if(__name__=='__main__'):
    main()    