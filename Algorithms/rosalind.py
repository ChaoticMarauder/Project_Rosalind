import numpy as np

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