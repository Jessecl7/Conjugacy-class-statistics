def conjstatistics(alpha, n):  #Takes in a conjugacy classes of permutations, then outputs various statistics about it.
    G = SymmetricGroup(int(n))
    alpha = G(alpha)
    Cclass = G.conjugacy_class(alpha)
    C = Cclass.list()
    clen = len(C)
    inv = [0]*(n^2)
    maj = [0]*(n^2)
    invprod = [0]*(n^2)
    majprod = [0]*(n^2)
    triples = 0
    pairs = 0
    for j in range(clen):
        temp = G(C[j])
        tempnum = 0
        temp = str(temp)
        temp = Permutation(temp)
        tempnum = temp.number_of_inversions()
        inv[tempnum] = inv[tempnum]+1
        tempnum = temp.major_index()
        maj[tempnum] = maj[tempnum]+1
    print("Average number of inversions in the first class")
    temp = 0
    for i in range(n^2):
        temp = temp + i*inv[i]
    show(temp/clen)
    print("Second moment")
    temp = 0
    for i in range(n^2):
        temp = temp + i^2*inv[i]
    show(temp/clen)
    print("Third moment")
    temp = 0
    for i in range(n^2):
        temp = temp + i^3*inv[i]
    show(temp/clen)
    print("Average major index of the first class")
    temp = 0
    for i in range(n^2):
        temp = temp + i*maj[i]
    show(temp/clen)
    print("Second moment")
    temp = 0
    for i in range(n^2):
        temp = temp + i^2*maj[i]
    show(temp/clen)
    print("Third moment")
    temp = 0
    for i in range(n^2):
        temp = temp + i^3*maj[i]
    show(temp/clen)
    temp = 0
    for i in range(n^2):
        temp = temp + i*invprod[i]
    temp = 0
    for i in range(n^2):
        temp = temp + i*majprod[i]
    print("Coefficients of the inversion polynomial are:")
    show(inv)
    print("Coefficients of the major index polynomial are:")
    show(maj)
    return