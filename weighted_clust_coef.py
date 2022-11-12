import numpy as np
import pandas as pd
import sys # to use input parameters and issue error

# Weighted clustering coefficient by Barrat et al. (2004)
def Cnet_barrat(nV, nE, E, w):

    k = np.zeros(nV) # degree
    s = np.zeros(nV) # node strength
    E = E.astype('int')
    for i in range(nE):
        k[E[i, 0]] += 1 # degree
        k[E[i, 1]] += 1
        s[E[i, 0]] += w[i] # node strength
        s[E[i, 1]] += w[i]

    for i in range(nE):
        if E[i, 0]==E[i, 1]:
            sys.exit("self-loop disallowed") # " << E[2*i] << endl;
        elif E[i, 0] > E[i, 1]:
            E[i, 0], E[i, 1] = E[i, 1], E[i, 0] # swap
        # now E[i, 0] < E[i, 1]

#    x = np.hstack((E, w))
#    np.savetxt('tmpx.txt', x, fmt=['%d', '%d', '%1.3f'])
    sorted_index = np.argsort(E[:,0])
    E = E[sorted_index]
    w = w[sorted_index]
#    y = np.hstack((E, w))
#    np.savetxt('tmpy.txt', y, fmt=['%d', '%d', '%1.3f'])

    P = np.zeros(nV, dtype=int) # P[i] = number of edges whose starting node's index <= i
    for i in range(nE):
            P[E[i, 0]] += 1
    for i in range(1, nV):
        P[i] = P[i] + P[i-1]

    localC = np.zeros(nV)
    # calculate the local weighted clustering coefficient
    for vs in range(nV):
        if vs==0:
            first = 0
        else:
            first = P[vs-1]

        for i in range(first, P[vs]):
            for j in range(i+1, P[vs]):
                ve1 = E[i, 1]; # ve1 is a neighbor of vs
                ve2 = E[j, 1]; # ve2 is another neighbor of vs. We have ve2 > ve1 guaranteed.
                if E[i, 0] != vs or E[j, 0] != vs:
                    sys.exit("vs inconsistent")
                if vs>=ve1 or vs>=ve2:
                    sys.exit("vs<ve1,ve2 violated")
                if ve1==ve2:
                    sys.exit("ve1 and ve2 must be different") #  << vs << " " << ve1 << endl;

                # examine if ve1 and ve2 are adjacent
                ii=P[ve1-1] # note that ve >= 1 is guaranteed
                found=0
                while ii<P[ve1] and found==0:
                    if E[ii, 1]==ve2:
                        found=1 # triangle
                    else:
                        ii += 1
                if found==0: # this part is necessary because np.argsort does not preserve the ascending order in E[:,1]
                    ii=P[ve2-1] # note that ve >= 1 is guaranteed
                    while ii<P[ve2] and found==0:
                        if E[ii, 1]==ve1:
                            found=1 # triangle
                        else:
                            ii += 1

                # edge i: (vs,ve1), edge j: (vs,ve2), edge ii: (ve1,ve2)
                if found==1:
                    localC[vs] += w[i] + w[j]
                    localC[ve1] += w[i] + w[ii]
                    localC[ve2] += w[j] + w[ii]
    # all triplets done

    localC = localC / (s * (k-1)) # Elementwise operation. Becomes inf if k[i]==0 and k[i]==1, but anyways
    C = np.average(localC[k>=2]) # global weighted clustering coefficient
    nV_eff = np.count_nonzero(k>=2) # number of nodes with degree >= 2
    if nV_eff < nV:
        print(nV - nV_eff, "nodes have degree 0 or 1")

    return C, localC

# Weighted clustering coefficient by Onnela et al. (2005)
def Cnet_onnela(nV, nE, E, w):

    k = np.zeros(nV) # degree
    s = np.zeros(nV) # node strength
    E = E.astype('int')
    for i in range(nE):
        k[E[i, 0]] += 1 # degree
        k[E[i, 1]] += 1
        s[E[i, 0]] += w[i] # node strength
        s[E[i, 1]] += w[i]

    for i in range(nE):
        if E[i, 0]==E[i, 1]:
            sys.exit("self-loop disallowed") # " << E[2*i] << endl;
        elif E[i, 0] > E[i, 1]:
            E[i, 0], E[i, 1] = E[i, 1], E[i, 0] # swap
        # now E[i, 0] < E[i, 1]

#    x = np.hstack((E, w))
#    np.savetxt('tmpx.txt', x, fmt=['%d', '%d', '%1.3f'])
    sorted_index = np.argsort(E[:,0])
    E = E[sorted_index]
    w = w[sorted_index]
#    y = np.hstack((E, w))
#    np.savetxt('tmpy.txt', y, fmt=['%d', '%d', '%1.3f'])

    P = np.zeros(nV, dtype=int) # P[i] = number of edges whose starting node's index <= i
    for i in range(nE):
            P[E[i, 0]] += 1
    for i in range(1, nV):
        P[i] = P[i] + P[i-1]

    localC = np.zeros(nV)
    # calculate the local weighted clustering coefficient
    for vs in range(nV):
        if vs==0:
            first = 0
        else:
            first = P[vs-1]

        for i in range(first, P[vs]):
            for j in range(i+1, P[vs]):
                ve1 = E[i, 1]; # ve1 is a neighbor of vs
                ve2 = E[j, 1]; # ve2 is another neighbor of vs. We have ve2 > ve1 guaranteed.
                if E[i, 0] != vs or E[j, 0] != vs:
                    sys.exit("vs inconsistent")
                if vs>=ve1 or vs>=ve2:
                    sys.exit("vs<ve1,ve2 violated")
                if ve1==ve2:
                    sys.exit("ve1 and ve2 must be different") #  << vs << " " << ve1 << endl;

                # examine if ve1 and ve2 are adjacent
                ii=P[ve1-1] # note that ve >= 1 is guaranteed
                found=0
                while ii<P[ve1] and found==0:
                    if E[ii, 1]==ve2:
                        found=1 # triangle
                    else:
                        ii += 1
                if found==0: # this part is necessary because np.argsort does not preserve the ascending order in E[:,1]
                    ii=P[ve2-1] # note that ve >= 1 is guaranteed
                    while ii<P[ve2] and found==0:
                        if E[ii, 1]==ve1:
                            found=1 # triangle
                        else:
                            ii += 1

                # edge i: (vs,ve1), edge j: (vs,ve2), edge ii: (ve1,ve2)
                if found==1:
                    tmp = (w[i] * w[j] * w[ii])**(1/3)
                    localC[vs] += tmp
                    localC[ve1] += tmp
                    localC[ve2] += tmp
    # all triplets done

    localC = localC / (k * (k-1) / 2 * np.amax(w)) # Elementwise operation. Becomes inf if k[i]==0 and k[i]==1, but anyways
    C = np.average(localC[k>=2]) # global weighted clustering coefficient
    nV_eff = np.count_nonzero(k>=2) # number of nodes with degree >= 2
    if nV_eff < nV:
        print(nV - nV_eff, "nodes have degree 0 or 1")

    return C, localC

# Weighted clustering coefficient by Zhang and Horvath (2005)
def Cnet_zhang(nV, nE, E, w):

    k = np.zeros(nV) # degree
    s = np.zeros(nV) # node strength
    E = E.astype('int')
    for i in range(nE):
        k[E[i, 0]] += 1 # degree
        k[E[i, 1]] += 1
        s[E[i, 0]] += w[i] # node strength
        s[E[i, 1]] += w[i]

    for i in range(nE):
        if E[i, 0]==E[i, 1]:
            sys.exit("self-loop disallowed") # " << E[2*i] << endl;
        elif E[i, 0] > E[i, 1]:
            E[i, 0], E[i, 1] = E[i, 1], E[i, 0] # swap
        # now E[i, 0] < E[i, 1]

#    x = np.hstack((E, w))
#    np.savetxt('tmpx.txt', x, fmt=['%d', '%d', '%1.3f'])
    sorted_index = np.argsort(E[:,0])
    E = E[sorted_index]
    w = w[sorted_index]
#    y = np.hstack((E, w))
#    np.savetxt('tmpy.txt', y, fmt=['%d', '%d', '%1.3f'])

    P = np.zeros(nV, dtype=int) # P[i] = number of edges whose starting node's index <= i
    for i in range(nE):
            P[E[i, 0]] += 1
    for i in range(1, nV):
        P[i] = P[i] + P[i-1]

    localC = np.zeros(nV)
    denom = np.zeros(nV) # denominator of the local clustering coefficient
    # calculate the local weighted clustering coefficient
    for vs in range(nV):
        if vs==0:
            first = 0
        else:
            first = P[vs-1]

        for i in range(first, P[vs]):
            for j in range(i+1, P[vs]):
                ve1 = E[i, 1]; # ve1 is a neighbor of vs
                ve2 = E[j, 1]; # ve2 is another neighbor of vs. We have ve2 > ve1 guaranteed.
                if E[i, 0] != vs or E[j, 0] != vs:
                    sys.exit("vs inconsistent")
                if vs>=ve1 or vs>=ve2:
                    sys.exit("vs<ve1,ve2 violated")
                if ve1==ve2:
                    sys.exit("ve1 and ve2 must be different") #  << vs << " " << ve1 << endl;

                # examine if ve1 and ve2 are adjacent
                ii=P[ve1-1] # note that ve >= 1 is guaranteed
                found=0
                while ii<P[ve1] and found==0:
                    if E[ii, 1]==ve2:
                        found=1 # triangle
                    else:
                        ii += 1
                if found==0: # this part is necessary because np.argsort does not preserve the ascending order in E[:,1]
                    ii=P[ve2-1] # note that ve >= 1 is guaranteed
                    while ii<P[ve2] and found==0:
                        if E[ii, 1]==ve1:
                            found=1 # triangle
                        else:
                            ii += 1

                # edge i: (vs,ve1), edge j: (vs,ve2), edge ii: (ve1,ve2)
                if found==1:
                    tmp = w[i] * w[j] * w[ii]
                    localC[vs] += tmp
                    localC[ve1] += tmp
                    localC[ve2] += tmp
                    denom[vs] += w[i] * w[j]
                    denom[ve1] += w[i] * w[ii]
                    denom[ve2] += w[j] * w[ii]
    # all triplets done

    localC = localC / (denom * np.amax(w)) # Elementwise operation. Becomes inf if k[i]==0 and k[i]==1, but anyways
    C = np.average(localC[k>=2]) # global weighted clustering coefficient
    nV_eff = np.count_nonzero(k>=2) # number of nodes with degree >= 2
    if nV_eff < nV:
        print(nV - nV_eff, "nodes have degree 0 or 1")

    return C, localC