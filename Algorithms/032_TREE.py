def connected_tree(n, edge_list):
    current_edges = len(edge_list)
    edges_needed = (n-1) - current_edges
    
    return edges_needed

def main():
    with open('datasets/rosalind_tree.txt') as input_file:
        input_data = input_file.read().strip().split('\n')
        
    n = int(input_data.pop(0))
    edge_list = list(map(int,edge.split()) for edge in input_data)
    
    edges_needed = connected_tree(n, edge_list)
    
    print(str(edges_needed))
    
    with open('solutions/rosalind_tree.txt', 'w') as output_file:
        output_file.write(str(edges_needed))
        
if(__name__=='__main__'):
    main()
    