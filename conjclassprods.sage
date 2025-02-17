import random
import numpy
from matplotlib import pyplot as plt

#This is a collection of various programs that calculate statistics on various sets of permutations.
#Most of the functions work with products of conjugacy classes of permutations.

def invproducts(j,k,n): #Calculates the distribution of number of inversions for a product of two permutatations with given numbers of descents.
    G = SymmetricGroup(int(n))
    elts = G.list()
    eltslen = len(elts)
    inv = [0]*(n^2)
    for i in range(eltslen):
        for l in range(eltslen):
            temp = G(elts[i])
            temp = str(temp)
            temp = Permutation(temp)
            temp2 = G(elts[l])
            temp2 = str(temp2)
            temp2 = Permutation(temp2)
            if temp.number_of_descents() == j and temp2.number_of_descents() == k:
                temp = G(temp)*G(temp2)
                temp = str(temp)
                temp = Permutation(temp)
                tempnum = temp.number_of_descents()
                inv[tempnum] = inv[tempnum]+1
    show(inv)
                
    
def conjprodstatistics(alpha, alphatwo, n): #Takes in two conjugacy classes, and outputs the distributions and average of various statistics on their product.
    G = SymmetricGroup(int(n))
    alpha = G(alpha)
    alphatwo = G(alphatwo)
    Cclass = G.conjugacy_class(alphatwo)
    C = Cclass.list()
    clen = len(C)
    Dclass = G.conjugacy_class(alpha)
    D = Dclass.list()
    dlen = len(D)
    inv = [0]*(n^2)
    maj = [0]*(n^2)
    invprod = [0]*(n^2)
    majprod = [0]*(n^2)
    for j in range(clen):
        temp = G(C[j])
        temp = str(temp)
        temp = Permutation(temp)
        tempnum = temp.number_of_inversions()
        inv[tempnum] = inv[tempnum]+1
        tempnum = temp.major_index()
        maj[tempnum] = maj[tempnum]+1
        for k in range(dlen):
            temp = G(D[k]) * G(C[j])
            temp = str(temp)
            temp = Permutation(temp)
            tempnum = temp.number_of_inversions()
            invprod[tempnum]=invprod[tempnum]+1
            tempnum = temp.major_index()
            majprod[tempnum]=majprod[tempnum]+1
    print("Average number of inversions in the first class")
    #show(inv)
    temp = 0
    for i in range(n^2):
        temp = temp + i*inv[i]
    show(temp/factorial(n-1))
    print("Average major index of the first class")
    #show(maj)
    temp = 0
    for i in range(n^2):
        temp = temp + i*maj[i]
    show(temp/factorial(n-1))
    print("Average number of inversions in the class product")
    #show(invprod)
    temp = 0
    for i in range(n^2):
        temp = temp + i*invprod[i]
    show(temp/(factorial(n-1))^2)
    print("Average major index of the class product")
    #show(majprod)
    temp = 0
    for i in range(n^2):
        temp = temp + i*majprod[i]
    show(temp/(factorial(n-1))^2)
    return
        
def expectedno(alpha, alphatwo, n): #Calculates the average number of cycles in a product of two conjugacy classes.
    G = SymmetricGroup(int(n))
    alpha = G(alpha)
    alphatwo = G(alphatwo)
    Cclass = G.conjugacy_class(alphatwo)
    C = Cclass.list()
    clen = len(C)
    D = G.conjugacy_class(alpha)
    alphazero = D.representative()
    totalcycles = 0
    for j in range(clen):
        temp = G(alphazero) * G(C[j])
        temp = temp.cycle_tuples(singletons=True)
        temptype = str(temp)
        tempnum = temptype.count("(")
        totalcycles = totalcycles + tempnum
    totalcycles = totalcycles/clen
    totalcyclesdec  = totalcycles.numerical_approx(digits=5)
    print(str(alphatwo) + " and " + str(alpha) + " have expected number of cycles " + str(totalcyclesdec))
          
    return

def conjugacyprod(alpha, alphatwo, n): #Calculates a the polynomial that counts the number of cycles in a conjugacy class product.  Performs a change of basis on this polynomial then calculates the roots. 
    G = SymmetricGroup(int(n))
    alpha = G(alpha)
    alphatwo = G(alphatwo)
    Cclass = G.conjugacy_class(alphatwo)
    C = Cclass.list()
    clen = len(C)
    D = G.conjugacy_class(alpha)
    alphazero = D.representative()
    totalcycles = 0
    dist = [0] * (n+1)
    for j in range(clen):
        temp = G(alphazero) * G(C[j])
        temp = temp.cycle_tuples(singletons=True)
        temp = str(temp)
        numcycles = temp.count('(')
        dist[numcycles] = dist[numcycles] + 1
    R.<t> = QQ[]
    poly = R(dist)
    show(poly)
    poly = poly - poly
    maxcoeff = numpy.max(numpy.nonzero(dist))
    eulers = [0] * (maxcoeff + 1)
    eulers[0] = 1
    for i in range(maxcoeff):
        i = i+1
        for k in range(i):
            eulers[i] = eulers[i] + binomial(i, k) * eulers[k] * (t-1)^(i - 1 - k)
    for i in range(maxcoeff):
        i = i+1
        poly = poly + dist[i] * (1-t)^(maxcoeff - i) * eulers[i]   
    poly = poly / factorial(int(maxcoeff))
    print(str(poly))
    show(poly.roots(CC))
    return

def cyclelengths(alpha, alphatwo, n): #Calculates the lengths of the different cycles, across all permutations in a product of conjugacy classes.
    G = SymmetricGroup(int(n))
    alpha = G(alpha)
    alphatwo = G(alphatwo)
    Cclass = G.conjugacy_class(alphatwo)
    C = Cclass.list()
    clen = len(C)
    D = G.conjugacy_class(alpha)
    alphazero = D.representative()
    totalcycles = 0
    lengths = [0] * n
    for j in range(clen):
        temp = G(alphazero) * G(C[j])
        temp = str(temp)
        temp = tuple(temp)
        if temp[1] == "1":
            for k in range(len(temp)):
                if temp[k] == ")":
                    lengths[k/2 - 1] = lengths[k/2-1] + 1
                    break
        else: lengths[0] = lengths[0] + 1
            
    show(lengths)
    return

         