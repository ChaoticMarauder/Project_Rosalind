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