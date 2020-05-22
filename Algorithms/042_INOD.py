def internal_nodes(n):
    # n is the number of leaves on the unrooted binary tree
    return n-2

def main():
    with open('datasets/rosalind_inod.txt') as input_file:
        n = int(input_file.read().strip())
    print(internal_nodes(n))
    
    with open('solutions/rosalind_inod.txt', 'w') as output_file:
        output_file.write(str(internal_nodes(n)))
if(__name__=='__main__'):
    main()