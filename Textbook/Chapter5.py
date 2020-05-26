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