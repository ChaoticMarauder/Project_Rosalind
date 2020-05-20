def string_slice(string, list_index):
    list_words=['','']
    list_words[0]=string[list_index[0]:list_index[1]+1]
    list_words[1]=string[list_index[2]:list_index[3]+1]
    
    return list_words

def main():
    
    with open('datasets/rosalind_ini3.txt') as input_file:
         string, list_index = [line.strip() for line in input_file.readlines()]
     
    list_index=list(map(int,list_index.split()))
     
    words=string_slice(string,list_index)
    print(' '.join(words))
     
    with open('solutions/rosalind_ini3.txt', 'w') as output_file:
        output_file.write(' '.join(words))
    
if(__name__=='__main__'):
    main()