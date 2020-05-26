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

def main():
    with open('datasets/rosalind_ba5a.txt') as input_file:
        money, coin_list = input_file.read().strip().split('\n')
        
    money = int(money)
    coin_list = coin_list.split(',')
    coin_list = list(map(int, coin_list))
    
   
    min_coins = dynamic_change(money, coin_list)
    
    print(str(min_coins))
    
    
    with open('solutions/rosalind_ba5a.txt', 'w') as output_file:
        output_file.write(str(min_coins)) 
        
if __name__ == '__main__':
    main()