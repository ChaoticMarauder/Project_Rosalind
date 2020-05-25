def linear_spectrum(peptide):
    dict_weights={'A':71,'C':103,'D':115,'E':129,'F':147,'G':57,'H':137,
                  'I':113,'K':128,'L':113,'M':131,'N':114,'P':97,'Q':128,
                  'R':156,'S':87,'T':101,'V':99,'W':186,'Y':163}
    
    prefix_mass = [0]
    
    for i in range(1,len(peptide)+1):
        for key in dict_weights:
            if(peptide[i-1] == key):
                prefix_mass.append(prefix_mass[i-1]+dict_weights[key])
    
    spectrum=[]            
    for i in range(len(peptide)):
        for j in range(i+1,len(peptide)+1):
            weight_sample = prefix_mass[j]-prefix_mass[i]
            spectrum.append(weight_sample)
            
    linear_spectrum_list = sorted(spectrum)
   
    return linear_spectrum_list 

def cyclic_spectrum(peptide):
    dict_weights={'A':71,'C':103,'D':115,'E':129,'F':147,'G':57,'H':137,
                  'I':113,'K':128,'L':113,'M':131,'N':114,'P':97,'Q':128,
                  'R':156,'S':87,'T':101,'V':99,'W':186,'Y':163}
    
    
    
    prefix_mass = [0]
    
    for i in range(1,len(peptide)+1):
        for key in dict_weights:
            if(peptide[i-1] == key):
                prefix_mass.append(prefix_mass[i-1]+dict_weights[key])
    
    
    peptide_mass = prefix_mass[len(peptide)]
    spectrum=[0]            
    for i in range(len(peptide)):
        for j in range(i+1,len(peptide)+1):
            weight_sample = prefix_mass[j]-prefix_mass[i]
            spectrum.append(weight_sample)
            if i>0 and j<len(peptide):
                weight_sample = peptide_mass-(prefix_mass[j]-prefix_mass[i])
                spectrum.append(weight_sample)
                
    cyclic_spectrum_list = sorted(spectrum)
    
    return cyclic_spectrum_list 

def main():
    with open('datasets/rosalind_ba4c.txt') as input_file:
        peptide = input_file.read().strip()
        
    cyclic_spectrum_list = cyclic_spectrum(peptide)
    
    print(' '.join(list(map(str,cyclic_spectrum_list))))
    
    with open('solutions/rosalind_ba4c.txt', 'w') as output_file:
        output_file.write(' '.join(list(map(str,cyclic_spectrum_list))))
        
if(__name__=='__main__'):
    main()