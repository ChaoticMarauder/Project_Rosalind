import numpy as np
def manhattan_tourist(n, m, Down, Right):
    s = np.zeros([n+1, m+1], dtype = int)
    
    for i in range(1,n+1):
        s[i][0] = s[i-1][0] + Down[i-1][0]
    
    for j in range(1,m+1):
        s[0][j] = s[0][j-1] + Right[0][j-1]
    
    for i in range(1,n+1):
        for j in range(1,m+1):
            
            s[i][j] = max(s[i-1][j]+Down[i-1][j],s[i][j-1]+Right[i][j-1])
    
            
    return s[n][m]

def main():
    with open('datasets/rosalind_ba5b.txt') as input_file:
        Down, Right = input_file.read().strip().split('-')
    
    Down = Down.strip().split('\n')
    Right = Right.strip().split('\n')    
    n, m = Down[0].split()
    Down = Down[1:]
    
    n = int(n)
    m = int(m)
    
    down_temp = []
    right_temp = []
    for i in range(n):
        a = Down[i].split()
        down_temp.append(list(map(int,a)))
    
    for j in range(n+1):
        a = Right[j].split()
        right_temp.append(list(map(int,a)))
    
    Down = down_temp
    Right = right_temp

   
    
    score = manhattan_tourist(n, m, Down, Right)
    
    print(str(score))
    
    
    with open('solutions/rosalind_ba5b.txt', 'w') as output_file:
        output_file.write(str(score)) 
        
if __name__ == '__main__':
    main()