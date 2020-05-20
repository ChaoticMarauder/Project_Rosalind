def dict_words(string):
    dict_word={}
    
    list_words=string.split()
    
    for word in list_words:
        if(dict_word.get(word)!=None):
            dict_word[word]+=1
        else:
            dict_word[word]=1
                
    
    return dict_word

def main():
    with open('datasets/rosalind_ini6.txt') as input_file:
        string=input_file.read().strip()
        
    dict_word=dict_words(string)
    
    for key, value in dict_word.items():
        print(str(key)+' '+str(value))
        
    with open('solutions/rosalind_ini6.txt', 'w') as output_file:
        for key, value in dict_word.items():
            output_file.write(str(key)+' '+str(value))
    
if(__name__=='__main__'):
    main()