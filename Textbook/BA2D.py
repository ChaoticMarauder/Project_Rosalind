import math

def SymbolToNumber(Symbol):

    if Symbol == "A":
        return 0
    elif Symbol == "C":
        return 1
    elif Symbol == "G":
        return 2
    elif Symbol == "T":
        return 3


def NumberToSymbol(index):

    if index == 0:
        return str("A")
    elif index == 1:
        return str("C")
    elif index == 2:
        return str("G")
    elif index == 3:
        return str("T")


def HammingDistance(p, q):

    return sum(s1 != s2 for s1, s2 in zip(p, q))


def window(s, k):
    for i in range(1 + len(s) - k):
        yield s[i:i+k]


def ProfileMostProbable(Text, k, Profile):

    letter = [[] for key in range(k)]
    probable = ""
    hamdict = {}
    index = 1
    for a in range(k):
        for j in "ACGT":
            letter[a].append(Profile[j][a])
    for b in range(len(letter)):
        number = max(letter[b])
        probable += str(NumberToSymbol(letter[b].index(number)))
    for c in window(Text, k):
        for x in range(len(c)):
            y = SymbolToNumber(c[x])
            index *= float(letter[x][y])
        hamdict[c] = index
        index = 1
    for pat, ham in hamdict.items():
        if ham == max(hamdict.values()):
            final = pat
            break
    return final


def Count(Motifs):

    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for i in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def FindConsensus(motifs):

    consensus = ""
    for i in range(len(motifs[0])):
        countA, countC, countG, countT = 0, 0, 0, 0
        for motif in motifs:
            if motif[i] == "A":
                countA += 1
            elif motif[i] == "C":
                countC += 1
            elif motif[i] == "G":
                countG += 1
            elif motif[i] == "T":
                countT += 1
        if countA >= max(countC, countG, countT):
            consensus += "A"
        elif countC >= max(countA, countG, countT):
            consensus += "C"
        elif countG >= max(countC, countA, countT):
            consensus += "G"
        elif countT >= max(countC, countG, countA):
            consensus += "T"
    return consensus


def ProfileMatrix(motifs):

    Profile = {}
    A, C, G, T = [], [], [], []
    for j in range(len(motifs[0])):
        countA, countC, countG, countT = 0, 0, 0, 0
        for motif in motifs:
            if motif[j] == "A":
                countA += 1
            elif motif[j] == "C":
                countC += 1
            elif motif[j] == "G":
                countG += 1
            elif motif[j] == "T":
                countT += 1
        A.append(countA)
        C.append(countC)
        G.append(countG)
        T.append(countT)
    Profile["A"] = A
    Profile["C"] = C
    Profile["G"] = G
    Profile["T"] = T
    return Profile


def Score(motifs):

    consensus = FindConsensus(motifs)
    score = 0.0000
    for motif in motifs:
        score += HammingDistance(consensus, motif)
    #print(score)
    return round(score, 4)

def GreedyMotifSearch(DNA, k, t):
    
    
    bestMotifs = []
    bestScore = math.inf
    for string in DNA:
        bestMotifs.append(string[:k])
    base = DNA[0]
    for i in window(base, k):
        # Change here. Should start with one element in motifs and build up.
        # As in the line "motifs ← list with only Dna[0](i,k)"
        # newMotifs = []
        newMotifs = [i]
        # Change here to iterate over len(DNA). 
        # Should go through "for j from 1 to |Dna| - 1"
        # for j in range(t):
        for j in range(1, len(DNA)):
            # Change here. Should build up motifs and build profile using them.
            # profile = ProfileMatrix([i])
            profile = ProfileMatrix(newMotifs)
            probable = ProfileMostProbable(DNA[j], k, profile)
            newMotifs.append(probable)

        # Change to < rather < = to ensure getting the most recent hit. As referenced in the instructions:
        # If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring **first**.
        if Score(newMotifs) < bestScore:
        #if Score(newMotifs) <= bestScore:
            bestScore = Score(newMotifs)
            bestMotifs = newMotifs
    return bestMotifs

def main():
    with open('datasets/rosalind_ba2d.txt') as input_data:
         A = input_data.read().strip().split('\n')
        
    k, t = A[0].split()
    
    k = int(k)
    t = int(t)
    
    DNA = A[1:]
    
    
    bestMotifs = GreedyMotifSearch(DNA, k, t)
    
    print('\n'.join(bestMotifs))
    
    with open('solutions/rosalind_ba2d.txt', 'w') as output_file:
        output_file.write('\n'.join(bestMotifs))
        
if(__name__=='__main__'):
    main()