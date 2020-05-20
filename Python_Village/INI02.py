def hypotenuse(a,b):
    return int(a**2+b**2)

def main():
     with open('datasets/rosalind_ini2.txt') as input_file:
         a,b = input_file.read().strip().split()
     
     h=hypotenuse(int(a),int(b))
     print(h)
     
     with open('solutions/rosalind_ini2.txt', 'w') as output_file:
         output_file.write(str(h))
    
if(__name__=='__main__'):
    main()