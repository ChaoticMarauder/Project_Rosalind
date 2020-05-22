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

def main():
    
    with open('datasets/rosalind_seto.txt') as input_file:
        U, A, B = input_file.readlines()
        
    U = set([i for i in range(1,int(U.strip())+1)])
    A = set(map(int, A.strip().rstrip('}').lstrip('{').replace(' ','').split(',')))
    B = set(map(int, B.strip().rstrip('}').lstrip('{').replace(' ','').split(',')))
    
    list_of_operations = set_operations(U, A, B)
    
    for idx in range(len(list_of_operations)):
        print('{'+', '.join(list(map(str,list_of_operations[idx])))+'}')
        
    with open('solutions/rosalind_seto.txt', 'w') as output_file:
        for idx in range(len(list_of_operations)):
            output_file.write('{'+', '.join(list(map(str,list_of_operations[idx])))+'}'+'\n')
            
if(__name__=='__main__'):
    main()
        