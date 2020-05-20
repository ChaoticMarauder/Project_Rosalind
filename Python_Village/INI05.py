def main():
    c=0
    line_list=[]
    with open('datasets/rosalind_ini5.txt') as input_file:
        for line in input_file.readlines():
            c=c+1
            if(c%2==0):
               line_list.append(line)
    
    print('\n'.join(line_list))
    
    with open('solutions/rosalind_ini5.txt', 'w') as output_file:
        output_file.write('\n'.join(line_list))
    
if(__name__=='__main__'):
    main()