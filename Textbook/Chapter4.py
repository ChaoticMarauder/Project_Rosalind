import math
from Chapter1 import reverse_complement

def translate(rna):
    dic = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
       "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
    protein=""
    l=len(rna)
    start_codon=0
    end_codon=1
    start=rna.find("AUG")
    if start!=-1:
        start_codon=start
        while start+2<l:
            i=rna[start:start+3]
            if (i=="UAA" or i=="UAG" or i=="UGA"):
                end_codon=start
                break
            for j in dic:
                if(j==i):
                    protein=protein+dic[i]
            start=start+3
    
    if(end_codon>1):
        return protein, start_codon, end_codon
    else:
        return "",start_codon, end_codon
    
    
def transcribe(dna):
    rna=[]
    for base in dna:
        if(base=='T'):
            rna.append('U')
        else:
            rna.append(base)
    
    return ''.join(rna)

def encode_peptide(dna, peptide):
    
    dict_codon = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"X", "UAG":"X",
       "UGU":"C", "UGC":"C", "UGA":"X", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
    
    rna = transcribe(dna)
    rev_comp = reverse_complement(dna)
    rev_comp_rna = transcribe(rev_comp)
    
    prot_candidates=['','','','','','']
    for i in range(int(len(rna)/3)):
        if(len(rna[3*i:3*i+3])==3):
            prot_candidates[0]=prot_candidates[0]+dict_codon[rna[3*i:3*i+3]]
        if(len(rna[3*i+1:3*i+4])==3):
            prot_candidates[1]=prot_candidates[1]+dict_codon[rna[3*i+1:3*i+4]]
        if(len(rna[3*i+2:3*i+5])==3):
            prot_candidates[2]=prot_candidates[2]+dict_codon[rna[3*i+2:3*i+5]]
        if(len(rev_comp_rna[3*i:3*i+3])==3):
            prot_candidates[3]=prot_candidates[3]+dict_codon[rev_comp_rna[3*i:3*i+3]]
        if(len(rev_comp_rna[3*i+1:3*i+4])==3):
            prot_candidates[4]=prot_candidates[4]+dict_codon[rev_comp_rna[3*i+1:3*i+4]]
        if(len(rev_comp_rna[3*i+2:3*i+5])==3):
            prot_candidates[5]=prot_candidates[5]+dict_codon[rev_comp_rna[3*i+2:3*i+5]]
        
    list_dna=[]    
    k = len(peptide)    
    for i in range(0,3):
        for j in range(len(prot_candidates[i])-k+1):
            if(prot_candidates[i][j:j+k] == peptide):
                pattern = dna[3*j+i:3*(j+k-1)+i+3]
                list_dna.append(pattern)
    
    for i in range(3,6):
        for j in range(len(prot_candidates[i])-k+1):
            if(prot_candidates[i][j:j+k] == peptide):
                pattern = rev_comp[3*j+i-3:3*(j+k-1)+i-3+3]
                list_dna.append(reverse_complement(pattern))       

    return list_dna

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

def counting_peptides(mass):
    dict_weights={'A':71,'C':103,'D':115,'E':129,'F':147,'G':57,'H':137,
                  'I':113,'K':128,'L':113,'M':131,'N':114,'P':97,'Q':128,
                  'R':156,'S':87,'T':101,'V':99,'W':186,'Y':163}
    
    list_weights=[]
    
    for key, value in dict_weights.items():
        if value not in list_weights:
            list_weights.append(value)
            
            
    len_list = len(list_weights) 
    
    count_peptides=[]
    
    for i in range(mass+1):
        count_peptides.append(0)
    
    count_peptides[0] = 1

    for i in range(1, mass+1):
        for j in range(len_list):
            if(i >= list_weights[j]):
                count_peptides[i] = count_peptides[i] + count_peptides[i-list_weights[j]]
                
    return count_peptides[mass]

def ProteinWeightDict():
	'''Returns a dictionary that translates Protein to Monoisotopic Mass.'''
	table ='''A   71.03711
	C   103.00919
	D   115.02694
	E   129.04259
	F   147.06841
	G   57.02146
	H   137.05891
	I   113.08406
	K   128.09496
	L   113.08406
	M   131.04049
	N   114.04293
	P   97.05276
	Q   128.05858
	R   156.10111
	S   87.03203
	T   101.04768
	V   99.06841
	W   186.07931
	Y   163.06333''' 

	protein_weight_dict = dict()

	for protein in table.split('\n'):
		protein_weight_dict[protein.strip('\t').split()[0]] = float(protein.strip('\t').split()[1])

	return protein_weight_dict

def append_char(add_list, add_chars):
	
	newlist = []
	for item in add_list:
		newlist += [item+ch for ch in set(add_chars)]
	return newlist

def spectrum(peptide):
	
	# Dictionary translating RNA to Protein
	weight = ProteinWeightDict()
	# Initialize as the mass 0 and the mass of the entire peptide.
	spec = [0, sum([int(weight[protein]) for protein in peptide])]
	# Find the masses of the adjacent intermediary subpeptides
	spec += [sum([int(weight[protein]) for protein in peptide[j:j+i]]) for i in range(1,len(peptide)) for j in range(len(peptide)-i+1)]
	# Sort the list in ascending order and convert to strings.
	spec = map(str,sorted(spec))

	return spec


def cyclo_peptide_sequencing(cyclospec):
    # Create the protein weight dictionary.
    weight = ProteinWeightDict()

    # Let n be the length of a given peptide, and L be the length of its cyclospectrum.  Then L = n(n-1) + 2.
    # Using the quadratic formula to to solve for n:  n = (sqrt(4L-7) + 1)/2
    n = int((math.sqrt(4*len(cyclospec)-7)+1)/2)
    protein, i = [], 1
    while(len(protein) != n):
        if int(cyclospec[i]) in map(int,weight.values()):
            protein.append(cyclospec[i])
        i+=1

    # Get the name of each protein corresponding to a given weight (if multiple, only take one).
    names = []
    for w in protein:
	    names.append([items[0] for items in weight.items() if int(items[1])==int(w)][0])

    # Build the possible sequences.
    seq = append_char(names,names)
    for repeat in range(1,n):
	    seq = filter(lambda subpeptide:set(spectrum(subpeptide)) < set(cyclospec), set(seq))
	    if repeat != n-1:
		    seq = append_char(seq,names)

    # Convert each protein to the proper format. 
    cyclopeptide_sequence = ['-'.join([str(int(weight[protein])) for protein in peptide]) for peptide in seq]
    
    return cyclopeptide_sequence