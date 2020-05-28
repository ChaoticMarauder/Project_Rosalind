from Bio.Seq import translate, CodonTable

with open('datasets/rosalind_ptra.txt') as input_data:
	coding_dna, protein = input_data.read().strip().split('\n')

for prot_id in CodonTable.ambiguous_generic_by_id.keys():
	if translate(coding_dna, table = prot_id, stop_symbol = '', to_stop=False) == protein:
		start_point = str(prot_id)
		break

print(start_point)
with open('solutions/rosalind_ptra.txt', 'w') as output_data:
	output_data.write(start_point)