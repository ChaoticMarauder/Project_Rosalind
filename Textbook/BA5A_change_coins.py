import math

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

def main():
    with open('datasets/rosalind_ba5a.txt') as input_file:
        money, coin_list = input_file.read().strip().split('\n')
        
    money = int(money)
    coin_list = coin_list.split(',')
    coin_list = list(map(int, coin_list))
    
   
    min_coins, change_coins = dynamic_change_coins(money, coin_list)
    
    print(str(min_coins))
    print(' '.join(list(map(str, change_coins))))
    print(str(sum(change_coins)))
    
    with open('solutions/rosalind_ba5a_change_coins.txt', 'w') as output_file:
        output_file.write(str(min_coins)+'\n') 
        output_file.write(' '.join(list(map(str, change_coins))))
if __name__ == '__main__':
    main()