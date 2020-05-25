import math

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

def main():
    with open('datasets/rosalind_ba4e.txt') as input_data:
        cyclospec = input_data.read().strip().split()
   
    
    cyclopeptide_sequence = cyclo_peptide_sequencing(cyclospec) 
    
    print(' '.join(cyclopeptide_sequence))
    
    with open('solutions/rosalind_ba4e.txt', 'w') as output_data:
	    output_data.write(' '.join(cyclopeptide_sequence))
        
if(__name__=='__main__'):
    main()
    
    