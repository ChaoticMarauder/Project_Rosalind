def gc(dna):
    length=len(dna)
    count=0
    for base in dna:
        if(base=='G' or base=='C'):
            count+=1
            
    gc_content=(count/length)*100
    
    return gc_content

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
                
       

def main():
    
    file_name='datasets/rosalind_gc.txt'
    
    seq_dict = parse_fasta(file_name)
    
    
    list_gc=[]
    list_id=[]
    
    for ID in seq_dict:
        dna_seq=seq_dict[ID]
        gc_content=gc(dna_seq)
        
        list_id.append(ID)
        list_gc.append(gc_content)
 
    init=0.0       
    max_index=0

    for index in range(len(list_gc)):
        if(list_gc[index]>=init):
            max_index=index
            init=list_gc[index]

    print(list_id[max_index])
    print(list_gc[max_index])
    
    with open('solutions/rosalind_gc.txt', 'w') as output_file:
        output_file.write(list_id[max_index]+'\n')
        output_file.write(str(list_gc[max_index]))
    

if(__name__=='__main__'):
    main()
    
            
        