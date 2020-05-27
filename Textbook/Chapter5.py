import math
import numpy as np

def dynamic_change(money, coin_list):
    minimum_coins = [0]
    
    for i in range(1, money+1):
        minimum_coins.append(math.inf)
        
        for j in range(len(coin_list)):
            if i>=coin_list[j]:
                
                if minimum_coins[i] >= minimum_coins[i-coin_list[j]]+1:
                    minimum_coins[i] = minimum_coins[i-coin_list[j]]+1
                    
    return minimum_coins[money]

def dynamic_change_coins(money, coin_list):
    minimum_coins = [0]
    coins = [0]
    for i in range(1, money+1):
        minimum_coins.append(math.inf)
        coins.append(0)
        for j in range(len(coin_list)):
            if i>=coin_list[j]:
                
                if minimum_coins[i] > minimum_coins[i-coin_list[j]]+1:
                    minimum_coins[i] = minimum_coins[i-coin_list[j]]+1
                    coins[i] = coin_list[j]
                    
    change_coins=[]
    value=money
    while(value!=0):
        change_coins.append(coins[value])
        value = value - coins[value]
    
    return minimum_coins[money], change_coins

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

def longest_common_subseq(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    
    s = np.zeros([m+1, n+1], dtype=int)
    
    for i in range(m):
        for j in range(n):
            if seq1[i] == seq2[j]:
                s[i+1][j+1] = s[i][j] + 1
            else:
                s[i+1][j+1] = max(s[i+1][j], s[i][j+1])
        
    lcs=''
    
    i, j = m, n
    
    while(i*j!=0):
        if s[i][j] == s[i-1][j]:
            i=i-1
        elif s[i][j] == s[i][j-1]:
            j=j-1    
        else:
            if seq1[i-1] == seq2[j-1]:
                lcs = seq1[i-1] + lcs
                i=i-1
                j=j-1
                
    return lcs

def global_alignment(seq1, seq2, scoring_matrix, indel_penalty):

    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]

    
    for i in range(1, m+1):
        s[i][0] = -i*indel_penalty
    for j in range(1, n+1):
        s[0][j] = -j*indel_penalty

    
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [s[i-1][j] - indel_penalty, s[i][j-1] - indel_penalty, s[i-1][j-1] + scoring_matrix[seq1[i-1], seq2[j-1]]]
            s[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(s[i][j])

    aligned_seq1 = seq1 
    aligned_seq2 = seq2

    
    max_score = s[m][n]
    i, j = m, n
    
    
    while(i*j != 0):
        if backtrack_matrix[i][j] == 0:
            aligned_seq2 = indel_inserted(aligned_seq2, j)
            i = i-1
        elif backtrack_matrix[i][j] == 1:
            aligned_seq1 = indel_inserted(aligned_seq1, i)
            j = j-1
        else:
            i = i-1
            j = j-1
            
    for indel in range(i):
        aligned_seq2 = indel_inserted(aligned_seq2, 0)
    for indel in range(j):
        aligned_seq1 = indel_inserted(aligned_seq1, 0)    

    return max_score, aligned_seq1, aligned_seq2

def indel_inserted(seq, i):
    return seq[:i] + '-' + seq[i:]

def local_alignment(seq1, seq2, scoring_matrix, indel_penalty):
    
    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    max_score = -10
    max_i = 0
    max_j = 0
    
        
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [0, s[i-1][j]-indel_penalty, s[i][j-1]-indel_penalty, s[i-1][j-1]+scoring_matrix[seq1[i-1], seq2[j-1]]]
            s[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(s[i][j])
            
            if s[i][j] > max_score:
                max_score = s[i][j]
                max_i = i
                max_j = j
                
    a = max_i
    b = max_j
    
    aligned_seq1 = seq1[:a]
    aligned_seq2 = seq2[:b]
    
    while(a*b!=0 and backtrack_matrix[a][b]!=0):
        if backtrack_matrix[a][b] == 1:
            aligned_seq2 = indel_inserted(aligned_seq2, b)
            a = a-1
        elif backtrack_matrix[a][b] == 2:
            aligned_seq1 = indel_inserted(aligned_seq1, a)
            b = b-1
        else:
            a = a-1
            b = b-1
            
    aligned_seq1 =  aligned_seq1[a:]
    aligned_seq2 =  aligned_seq2[b:]
    
    return max_score, aligned_seq1, aligned_seq2

def edit_distance(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    
    edit_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    for i in range(1, m+1):
        edit_matrix[i][0] = i
        
    for j in range(1, n+1):
        edit_matrix[0][j] = j
        
    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:
                edit_matrix[i][j] = edit_matrix[i-1][j-1]
            else:
                edit_matrix[i][j] = min(edit_matrix[i-1][j]+1, edit_matrix[i][j-1]+1, edit_matrix[i-1][j-1]+1)
    
    
    edit_distance_value = edit_matrix[m][n]
    
    return edit_distance_value

def fitting_alignment(seq1, seq2):
    
    m = len(seq1)
    n = len(seq2)
    
    f = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [f[i-1][j]-1, f[i][j-1]-1, f[i-1][j-1] + [-1,1][seq1[i-1]==seq2[j-1]]]
            f[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(f[i][j])
            
    
    i = max(enumerate([f[k][j] for k in range(n, m)]),key=lambda x: x[1])[0] + n
    j = n
    
    aligned_seq1 = seq1[:i]
    aligned_seq2 = seq2[:j]
    
    max_score = f[i][j]
    
    while(i*j != 0):
        if backtrack_matrix[i][j] == 0:
            i = i-1
            aligned_seq2 = indel_inserted(aligned_seq2, j)
        elif backtrack_matrix[i][j] == 1:
            j = j-1
            aligned_seq1 = indel_inserted(aligned_seq1, i)
        else:
            i=i-1
            j=j-1
            
    aligned_seq1 = aligned_seq1[i:]
    
    return max_score, aligned_seq1, aligned_seq2

def overlap_alignment(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    
    s = [[0 for i in range(n+1)] for j in range(m+1)]
    backtrack_matrix = [[0 for i in range(n+1)] for j in range(m+1)]
    
    max_score = -2*(m+n)
    max_i = 0
    max_j = 0
    
        
    for i in range(1, m+1):
        for j in range(1, n+1):
            score = [s[i-1][j]-2, s[i][j-1]-2, s[i-1][j-1] + [-2,1][seq1[i-1]==seq2[j-1]]]
            s[i][j] = max(score)
            backtrack_matrix[i][j] = score.index(s[i][j])
            
            if i == m or j == n:
                if s[i][j] > max_score:
                    max_score = s[i][j]
                    max_i = i
                    max_j = j
                
    a = max_i
    b = max_j
    
    aligned_seq1 = seq1[:a]
    aligned_seq2 = seq2[:b]
    
    while(a*b!=0):
        if backtrack_matrix[a][b] == 0:
            aligned_seq2 = indel_inserted(aligned_seq2, b)
            a = a-1
        elif backtrack_matrix[a][b] == 1:
            aligned_seq1 = indel_inserted(aligned_seq1, a)
            b = b-1
        else:
            a = a-1
            b = b-1
            
    aligned_seq1 =  aligned_seq1[a:]
    aligned_seq2 =  aligned_seq2[b:]
    
    return max_score, aligned_seq1, aligned_seq2