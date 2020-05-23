def SymbolToNumber(base):
    if(base=='A'):
        return 0
    if(base=='C'):
        return 1
    if(base=='G'):
        return 2
    if(base=='T'):
        return 3
    
def PatternToNumber(pattern):
    if(len(pattern)==0):
        return 0
    
    last_symbol = pattern[-1]
    prefix = pattern[0:len(pattern)-1]
    
    return 4*PatternToNumber(prefix)+SymbolToNumber(last_symbol)

def main():
    with open('datasets/rosalind_ba1l.txt') as input_file:
        pattern = input_file.read().strip()
        
    index = PatternToNumber(pattern)
    
    print(str(index))
    
    with open('solutions/rosalind_ba1l.txt', 'w') as output_file:
        output_file.write(str(index))

if(__name__=='__main__'):
    main()
    