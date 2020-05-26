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