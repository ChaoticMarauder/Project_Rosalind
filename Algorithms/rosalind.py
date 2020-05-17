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
    start=rna.find("AUG")
    if start!=-1:
        while start+2<l:
            i=rna[start:start+3]
            if (i=="UAA" or i=="UAG" or i=="UGA"):
                break
            for j in dic:
                if(j==i):
                    protein=protein+dic[i]
            start=start+3
    return protein


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