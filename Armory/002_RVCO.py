from Bio import SeqIO

with open('datasets/rosalind_rvco.txt') as input_data:
	rev_comp_matchings = sum([str(dna.seq) == str(dna.reverse_complement().seq) for dna in SeqIO.parse(input_data, 'fasta')])
	
print(str(rev_comp_matchings))
	
with open('solutions/rosalind_rvco.txt', 'w') as output_data:
    output_data.write(str(rev_comp_matchings))