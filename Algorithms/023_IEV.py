def expected_offspring(list_couples, num_offspring):
    expected_dominant = (num_offspring*list_couples[0] + num_offspring*list_couples[1] 
                         + num_offspring*list_couples[2] + 0.75*num_offspring*list_couples[3]
                         + 0.5*num_offspring*list_couples[4])
    
    return expected_dominant

def main():
    with open('datasets/rosalind_iev.txt') as input_file:
        couple_list = list(map(int,input_file.read().strip().split()))
        
    num=2    
    expected_dominant=expected_offspring(couple_list, num)
    
    print(expected_dominant)
    
    with open('solutions/rosalind_iev.txt', 'w') as output_file:
        output_file.write(str(expected_dominant))
        
if(__name__=='__main__'):
    main()