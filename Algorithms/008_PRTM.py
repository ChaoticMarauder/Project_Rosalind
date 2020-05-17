def protein_weight(protein):
    dict_weights={'A':71.03711,'C':103.00919,'D':115.02694,'E':129.04259,
              'F':147.06841,'G':57.02146,'H':137.05891,'I':113.08406,
              'K':128.09496,'L':113.08406,'M':131.04049,'N':114.04293,
              'P':97.05276,'Q':128.05858,'R':156.10111,'S':87.03203,
              'T':101.04768,'V':99.06841,'W':186.07931,'Y':163.06333}
    
    weight=0.0
    
    for amino_acid in protein:
        weight=weight+dict_weights[amino_acid]
        
    return weight   

def main():
    with open('datasets/rosalind_prtm.txt') as input_data:
        protein=input_data.read().strip()
    
    weight=round(protein_weight(protein),3)
    
    print(weight)
    
    with open('solutions/rosalind_prtm.txt', 'w') as output_file:
        output_file.write(str(weight))
        
    
if(__name__=='__main__'):
    main()