import numpy as np
from itertools import permutations
from itertools import product
from scipy.special import comb
import math

#parsing a fasta file
def parse_fasta(input_file):
    ID = None
    sequences = dict()
    
    with open(input_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            
            if (line[0]=='>'):
                ID = line[1:]
                sequences[ID] = ''
            else:
                sequences[ID] = sequences[ID] + line
                
    #returns a dictionary of sequence IDs as keys and sequences as values 
    return sequences

def gc(dna):
    length=len(dna)
    count=0
    for base in dna:
        if(base=='G' or base=='C'):
            count+=1
            
    gc_content=(count/length)*100
    
    return gc_content

def transcribe(dna):
    rna=[]
    for base in dna:
        if(base=='T'):
            rna.append('U')
        else:
            rna.append(base)
    
    return rna

def dna_base_count(dna):
    count_list=[]
    count_A=0
    count_C=0
    count_G=0
    count_T=0
    for base in dna:
        if(base=='A'):
            count_A+=1
        if(base=='C'):
            count_C+=1
        if(base=='G'):
            count_G+=1
        if(base=='T'):
            count_T+=1
    
    count_list=[count_A, count_C, count_G, count_T]
    
    return count_list

def reverse_complement(dna):
    rc=[]
    for base in dna:
        if(base=='T'):
            rc.append('A')
        if(base=='A'):
            rc.append('T')
        if(base=='G'):
            rc.append('C')
        if(base=='C'):
            rc.append('G')    
    
    
    rc.reverse()
    return rc

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


def hamming_distance(s,t):
    count=0
    if(len(s)==len(t)):
        for i in range(len(s)):
            if(s[i]!=t[i]):
                count+=1
    else:
        return -1
        
    return count

def dna_motif(dna, motif):
    m_len=len(motif)
    list_motifs=[]
    for i in range(len(dna)-m_len):
        if(dna[i:i+m_len]==motif):
            list_motifs.append(i+1)
            
    return list_motifs

def mortalfibonacciRabbits(n, m):
    rabbits = [0, 1, 1]
    for i in range(3, n + 1):
        if(i <= m):
            total = rabbits[i - 1] + rabbits[i - 2]
        elif(i == m + 1):
            total = rabbits[i - 1] + rabbits[i - 2] - 1
        else:
            total = rabbits[i - 1] + rabbits[i - 2] - rabbits[i - m - 1]
        rabbits.append(total)
    return (rabbits[n])

def rabbit_pairs(n,k):
    if(n==1):
        return 1
    if(n==2):
        return 1
    else:
        return k*rabbit_pairs(n-2,k)+rabbit_pairs(n-1,k)
    
    
def consensus_profile(seq_list):
    seq_length=len(seq_list[0])
    consensus_array=np.zeros([4,seq_length])
    
    for i in range(len(seq_list)):
        for j in range(seq_length):
            if(seq_list[i][j]=='A'):
                consensus_array[0][j]+=1
            if(seq_list[i][j]=='C'):
                consensus_array[1][j]+=1
            if(seq_list[i][j]=='G'):
                consensus_array[2][j]+=1
            if(seq_list[i][j]=='T'):
                consensus_array[3][j]+=1
                
    consensus_array=consensus_array.astype(int)
                
    dict_consensus={'A':consensus_array[0], 'C':consensus_array[1],
                        'G':consensus_array[2], 'T':consensus_array[3] }
        
    return dict_consensus, consensus_array


def base_mapping(n):
    if(n==0):
        base='A'
    if(n==1):
        base='C'
    if(n==2):
        base='G'
    if(n==3):
        base='T'
        
    return base
    
def consensus_sequence(array_cons):
    
    consensus_seq=""
    
    for i in range(np.size(array_cons,1)):
        max_value=0
        index=0
        for j in range(np.size(array_cons,0)):
            if(array_cons[j][i]>max_value):
                max_value=array_cons[j][i]
                index=j
            
        consensus_seq=consensus_seq+base_mapping(index)
            
    return consensus_seq

def reverse_palindrome(dna):
    if(dna==''.join(reverse_complement(dna))):
        return 1
    else:
        return 0
    
def restriction_sites(dna):
    len_dna=len(dna)
    
    site_list=[]
    length_list=[]
    
    for j in range(4,13):
        for i in range(len_dna-j+1):
                if(reverse_palindrome(dna[i:i+j])==1):
                    site_list.append(i+1)
                    length_list.append(j)
                    
    return site_list, length_list


def last_stop_codon(rna):
    
    last_stop=0
    for i in range(0,len(rna)-3,3):
        codon=rna[i:i+3]
        if(codon=="UAA" or codon=="UAG" or codon=="UGA"):
            last_stop=i+2
            
    return last_stop
        
def orf_translate(dna):
    
    rna=''.join(transcribe(dna))
    reverse_rna=''.join(transcribe(''.join(reverse_complement(dna))))
    
    orf_list=[]
    
    last_stop_rna=[]
    last_stop_reverse=[]
    
    for i in range(3):
        last_stop_rna.append(last_stop_codon(rna[i:]))
        last_stop_reverse.append(last_stop_codon(reverse_rna[i:]))
    
    
    first_strand=[]
    second_strand=[]
    
    
    for i in range(3):
        first_strand.append(rna[i:last_stop_rna[i]+i+1])
        second_strand.append(reverse_rna[i:last_stop_reverse[i]+i+1])
        
   
    
    for i in range(3):
        for j in range(len(first_strand[i])):
            seq1, start, end = translate(first_strand[i][j:])
            if(len(seq1)!=0):
                orf_list.append(seq1)

                
        for j in range(len(second_strand[i])):
            seq2, start, end = translate(second_strand[i][j:])
            if(len(seq2)!=0):
                orf_list.append(seq2)
            
    possible_orfs = [] 
    
    for seq in orf_list: 
        if seq not in possible_orfs:
            possible_orfs.append(seq) 
            
            
    return possible_orfs


def rna_splicing(dna, intron_list):
    
    for intron in intron_list:
        dna=dna.replace(intron,'')
    
    spliced_rna = ''.join(transcribe(dna))
    
    return spliced_rna


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
    
#genome assembly by generating the longest super string from given reads    
def max_overlap_sequence(seq_list):
    
    max_overlap=-1
    length_list=len(seq_list)
    
    for first_id in range(length_list):
        for second_id in range(length_list):
            if(first_id==second_id):
                continue
            else:
                first_string, second_string = seq_list[first_id], seq_list[second_id] 
                
                k=0
                
                while(first_string[k:]!=second_string[0:len(first_string[k:])]):
                    k=k+1
                
                if(len(first_string)-k > max_overlap):
                    max_overlap=len(first_string)-k
                    max_index=[first_id, second_id]
                    
    return [seq_list[max_index[0]]+seq_list[max_index[1]][max_overlap:]]+[seq_list[i] for i in range(length_list) if i not in max_index]

def shortest_super_sequence(seq_list):
    
    while(len(seq_list)>1):
        seq_list=max_overlap_sequence(seq_list)
        
    return seq_list[0]
#end of genome assembly

def factorial(n): 
    if(n<=1):
        return 1
    else:
        return n*factorial(n-1)
    
def rearrangements(n):
    num_set=range(1,n+1)
    permutation_list = list(permutations(num_set))
    
    return permutation_list
    
def k_mer_lexicographic(bases, k):
    k_mers = [''.join(item) for item in product(bases, repeat=k)]
    return k_mers

def signed_permutations(n):
    signed_ints_list=[]
    
    for i in range(1,n+1):
        signed_ints_list.append([i,-1*i])
    
    product_list=map(list,list(product(*signed_ints_list)))
    
    signed_permutations_list=[]
    for pro in product_list:
        signed_permutations_list+=map(list,list(permutations(pro)))
        
    return signed_permutations_list

def k_mer_sequence(dna):
    kmer_list=k_mer_lexicographic('ACGT',4)
    kmer_count=np.zeros(len(kmer_list))
    
    for i in range(len(dna)-3):
        kmer=dna[i:i+4]
        kmer_count[kmer_list.index(kmer)]+=1
    
    kmer_count=kmer_count.astype(int)
    return kmer_count

def mendel_first_law(dominant, heterozygous, recessive):
    num_offspring = comb(dominant+heterozygous+recessive, 2)
    
    num_recessive_offspring = (comb(recessive, 2) + 0.5*recessive*heterozygous
                              + 0.25*comb(heterozygous, 2))
    
    num_dominant_offspring = (comb(dominant, 2) + 0.75*comb(heterozygous, 2) 
                             + dominant*heterozygous + dominant*recessive 
                             + 0.5*recessive*heterozygous)
    
    ratio_recessive = num_recessive_offspring/num_offspring
    
    ratio_dominant = 1-ratio_recessive
    ratio_dominant = num_dominant_offspring/num_offspring
    
    return ratio_dominant

def expected_offspring(list_couples, num_offspring):
    expected_dominant = (num_offspring*list_couples[0] + num_offspring*list_couples[1] 
                         + num_offspring*list_couples[2] + 0.75*num_offspring*list_couples[3]
                         + 0.5*num_offspring*list_couples[4])
    
    return expected_dominant

def independent_alleles(k, N):
    prob_N=0.0
    total_offspring=2**k
    
    for i in range(N, total_offspring):
        #calculating the success of at least N in among 2**k trials-->binomial distribution
        prob_N=prob_N+comb(total_offspring,i)*((0.25)**i)*((0.75)**(total_offspring-i))
    
    prob_N=round(prob_N,4)
    return prob_N

def lex_variable(list_letters, n):
    final_list_letters=['$'] + list_letters
    
    kmer_list=list(product(final_list_letters, repeat=n))
    kmer_lex_list=[]
    for kmer in kmer_list:
        
        if('$' not in kmer):
            kmer_lex_list.append(''.join(kmer))
        
        else:
            for i in range(1,n):
                if(''.join(kmer[i:n])=='$'*(n-i) and '$' not in kmer[:i]):
                    kmer_lex_list.append(''.join(kmer).replace('$',''))
    
    return kmer_lex_list

def check_motif(motif, seq_list):
    for seq in seq_list:
        if((motif not in seq) or (len(motif)>len(seq))):
            return False
    
    return True
        
def longest_shared_motif(seq_list):
    
    longest_motif=''
    length_first=len(seq_list[0])
    longest_motif_length=0
    
    for i in range(length_first):
        for j in range(length_first, i, -1):
            if(j-i<longest_motif_length):
                break
            elif(check_motif(seq_list[0][i:j],seq_list)):
                longest_motif=seq_list[0][i:j]
                longest_motif_length=j-i
    
    return longest_motif


def prob_gc_seq(seq, gc_content):
    prob=0.0
    for i in range(len(seq)):
        
        if(seq[i]=='G' or seq[i]=='C'):
            prob=prob+math.log10((gc_content/2))
        
        if(seq[i]=='A' or seq[i]=='T'):
            prob=prob+math.log10((1-gc_content)/2)
            
    log_prob=round(prob,4)
    
    return log_prob

def list_prob_gc_seq(seq, list_gc_content):
    list_log_prob=[]
    for i in range(len(list_gc_content)):
        gc_content=list_gc_content[i]
        list_log_prob.append(prob_gc_seq(seq,gc_content))
        
    return list_log_prob

def matching_random_motif(N, gc_content, seq):
    GC_count=0
    AT_count=0
    
    for i in range(len(seq)):
        if(seq[i]=='G' or seq[i]=='C'):
            GC_count+=1
        else:
            AT_count+=1
            
    prob_seq = (gc_content/2)**GC_count*((1-gc_content)/2)**AT_count
    prob_not_motif = (1-prob_seq)**N
    
    prob_motif = 1-prob_not_motif
    
    return prob_motif

def motif_prob(seq, gc_content):
    GC_count=0
    AT_count=0
    
    for i in range(len(seq)):
        if(seq[i]=='G' or seq[i]=='C'):
            GC_count+=1
        else:
            AT_count+=1
            
    prob_seq = (gc_content/2)**GC_count*((1-gc_content)/2)**AT_count
    
    return prob_seq
    
def expected_restriction_sites(N, seq, gc_content):
    #seq is the restriction site sequence
    n = N - len(seq) + 1
    motif_probability = motif_prob(seq, gc_content)
    expected_sites = n*motif_probability
    
    return expected_sites

def mrna_to_protein(protein):
    codon_dict = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
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
    
    mrna_possibilities = list(codon_dict.values()).count("STOP")

    for aa in protein:
        mrna_possibilities = (mrna_possibilities*list(codon_dict.values()).count(aa))% 1000000

    return mrna_possibilities

def partial_perm(n, k):
    num_perm=1
    for i in range(n, n-k, -1):
        num_perm = (num_perm*i)%1000000
        
    return num_perm

def connected_tree(n, edge_list):
    current_edges = len(edge_list)
    edges_needed = (n-1) - current_edges
    
    return edges_needed

def distance_matrix(seq_list):
    num_seq=len(seq_list)
    length_seq=len(seq_list[0])
    distance_plot=np.zeros([num_seq, num_seq])
    
    for i in range(num_seq):
        for j in range(num_seq):
            distance_plot[i][j] = hamming_distance(seq_list[i], seq_list[j])/length_seq
            
    return distance_plot

def num_subsets(n):
    sset=1
    for i in range(n):
        sset=(sset*2)%1000000
    
    return sset

def alternate_splice(n,m):
    alternate_splice_count=0
    for k in range(m, n+1):
        alternate_splice_count = (alternate_splice_count + comb(n, k, exact=True))%1000000
        
    return alternate_splice_count

def independent_segregation(n):
    
    probability = 2**(2*n)
    
    log_prob_array = [-2*n*math.log10(2)]*2*n
    
    for k in range(2*n):
        probability = probability - comb(2*n, k, exact=True)
        log_prob_array[k]=log_prob_array[k] + math.log10(probability)
        
    return log_prob_array

def hardy_weinberg(p_list):
    recessive_list=[]
    
    for idx in range(len(p_list)):
        r = math.sqrt(p_list[idx])
        h = 2-r
        
        recessive_allele_probability = r*h
        recessive_list.append(recessive_allele_probability)
        
    return recessive_list

def set_operations(U, A, B):
    list_of_operations=[]
    
    union = list(A|B)
    list_of_operations.append(union)
    
    intersection = list(A&B)
    list_of_operations.append(intersection)
    
    difference_1 = list(A-B)
    list_of_operations.append(difference_1)
    
    difference_2 = list(B-A)
    list_of_operations.append(difference_2)
    
    A_complement = list(U-A)
    list_of_operations.append(A_complement)
    
    B_complement = list(U-B)
    list_of_operations.append(B_complement)
    
    return list_of_operations

def sex_linked_inheritance(allele_frequency_list):
    carrier_list=[]
    for idx in range(len(allele_frequency_list)):
        carrier = round(2*allele_frequency_list[idx]*(1-allele_frequency_list[idx]),3)
        carrier_list.append(carrier)
    
    return carrier_list

def de_bruijn_graph(seq_list):
    
    final_seq_list=[]
    
    for seq in seq_list:
        if seq not in final_seq_list:
            final_seq_list.append(seq)
    
    reverse_complement_seq_list=[]
    
    for seq in final_seq_list:
        reverse_complement_seq_list.append(''.join(reverse_complement(seq)))
        
    adjacency_list=[]
    k=len(final_seq_list[0])
    
    for seq in final_seq_list:
        edge=['','']
        edge[0]=seq[0:k-1]
        edge[1]=seq[1:k]
        if edge not in adjacency_list:
            adjacency_list.append(edge)
    
    for seq in reverse_complement_seq_list:
        edge=['','']
        edge[0]=seq[0:k-1]
        edge[1]=seq[1:k]
        if edge not in adjacency_list:
            adjacency_list.append(edge)
            
    return adjacency_list