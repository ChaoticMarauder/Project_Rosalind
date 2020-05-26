import math

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

