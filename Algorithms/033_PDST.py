from rosalind import parse_fasta
from rosalind import hamming_distance
import numpy as np

def distance_matrix(seq_list):
    num_seq=len(seq_list)
    length_seq=len(seq_list[0])
    distance_plot=np.zeros([num_seq, num_seq])
    
    for i in range(num_seq):
        for j in range(num_seq):
            distance_plot[i][j] = hamming_distance(seq_list[i], seq_list[j])/length_seq
            
    return distance_plot

def main():
    file_name='datasets/rosalind_pdst.txt'
    seq_dict = parse_fasta(file_name)
    
    seq_list=[]
    for key in seq_dict:
        seq_list.append(seq_dict[key])
        
    distance_plot=distance_matrix(seq_list)
    
    for i in range(len(distance_plot[0])):
        for j in range(len(distance_plot[0])):
            print(str(distance_plot[i][j])+' ')
        print()
    
    with open('solutions/rosalind_pdst.txt', 'w') as output_file:
        for i in range(len(distance_plot[0])):
            for j in range(len(distance_plot[0])):
                output_file.write(str(distance_plot[i][j])+' ')
            output_file.write('\n')
    
    
if(__name__=='__main__'):
    main()
            