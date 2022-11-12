import numpy as np
import pandas as pd
import math # to use constant pi
import sys # to use input parameters and issue error

# clustering coef for correlation matrices based on the partial correlation coefficient
# Input
#     cov: covariance matrix, whose (i,j) element is stored in cov[nV*i+j]
# Output
#     Ci_{type}: local clustering coefficient. An array of size nV should be supplied for storing the results.
#     C_{type}: global clustering coefficient. The pointer to a variable should be supplied for storing the results.
#     type = P, H, or A
#     P: C^{cor,P}. "P" for positive. Consider only positive edges.
#     H: C^{cor,H}. "H" for Heider. Consider all triangles and use the Heider's balance theory criterion.
#     A: C^{cor,A}. "A" for absolute value. Use the absolute value of the partial correlation coefficient.
#
def corrmat_C_PHA(cov):
    nV = len(cov) # number of nodes
    print(nV, 'nodes')
    corr = np.zeros((nV, nV))

    # calculate correlation matrix from covariance matrix
    for vs in range (nV):
        corr[vs, vs] = 1.0 # unused though
        for ve1 in range(vs):
            corr[vs, ve1] = cov[vs, ve1] / np.sqrt(cov[vs, vs] * cov[ve1, ve1])
            corr[ve1, vs] = corr[vs, ve1]

# The following is more concise but probably slower.
#
#    D_to_neg_half = np.diag(np.diag(cov) ** -0.5)
#    corr = np.dot(np.dot(D_to_neg_half, cov), D_to_neg_half)

    # initialization
    Ci_P = np.zeros(nV)
    Ci_H = np.zeros(nV)
    Ci_A = np.zeros(nV)
    denom_P = np.zeros(nV)
    denom_HA = np.zeros(nV)

    # local clustering coef
    for vs in range(nV):
        for ve1 in range(nV):
            if ve1 != vs:
                for ve2 in range(ve1):
                    if ve2 != vs:
                        # partial corr between ve1 and ve2 controlling for vs
                        pcorr = (corr[ve1, ve2] - corr[vs, ve1] * corr[vs, ve2]) / np.sqrt( (1 - corr[vs, ve1] * corr[vs, ve1]) * (1 - corr[vs, ve2] * corr[vs, ve2]) )
                        tmp = corr[vs, ve1] * corr[vs, ve2]
                        Ci_H[vs] += tmp * pcorr # Heider
                        Ci_A[vs] += np.absolute(tmp * pcorr)
                        denom_HA[vs] += np.absolute(tmp)
                        if corr[vs, ve1] > 0 and corr[vs, ve2] > 0: # (vs,ve1) and (vs,ve2) are both positive edges
                            Ci_P[vs] += tmp * pcorr
                            denom_P[vs] += tmp

    Ci_P /= denom_P
    Ci_H /= denom_HA
    Ci_A /= denom_HA
    C_P = np.average(Ci_P)
    C_H = np.average(Ci_H)
    C_A = np.average(Ci_A)

    return C_P, C_H, C_A, Ci_P, Ci_H, Ci_A

# Clustering coef for correlation matrices, C^{cor,A}, but only using +++ triangles or --- triangles
#
# A_pos: +++ triangles only
# A_neg: --- triangles only
# A: for absolute value. Use the absolute value of the partial correlation coefficient
def C_corA_pmonly(cov):
    nV = len(cov) # number of nodes
    print(nV, 'nodes')
    corr = np.zeros((nV, nV))

    # calculate correlation matrix from covariance matrix
    for vs in range (nV):
        corr[vs, vs] = 1.0 # unused though
        for ve1 in range(vs):
            corr[vs, ve1] = cov[vs, ve1] / np.sqrt(cov[vs, vs] * cov[ve1, ve1])
            corr[ve1, vs] = corr[vs, ve1]

# The following is more concise but probably slower.
#
#    D_to_neg_half = np.diag(np.diag(cov) ** -0.5)
#    corr = np.dot(np.dot(D_to_neg_half, cov), D_to_neg_half)

    # initialization
    localC_A_pos = np.zeros(nV)
    localC_A_neg = np.zeros(nV)
    denom_A_pos = np.zeros(nV)
    denom_A_neg = np.zeros(nV)

    # local clustering coef
    for vs in range(nV):
        for ve1 in range(nV):
            if ve1 != vs:
                for ve2 in range(ve1):
                    if ve2 != vs:
                        # partial corr between ve1 and ve2 controlling for vs
                        pcorr = (corr[ve1, ve2] - corr[vs, ve1] * corr[vs, ve2]) / np.sqrt( (1 - corr[vs, ve1] * corr[vs, ve1]) * (1 - corr[vs, ve2] * corr[vs, ve2]) )
                        if corr[vs, ve1] > 0 and corr[vs, ve2] > 0 and corr[ve1, ve2] > 0: # +++ triangle
                            tmp = corr[vs, ve1] * corr[vs, ve2] # This is positive
                            localC_A_pos[vs] += tmp * np.absolute(pcorr)
                            denom_A_pos[vs] += tmp
                        elif corr[vs, ve1] < 0 and corr[vs, ve2] < 0 and corr[ve1, ve2] < 0: # --- triangle
                            tmp = corr[vs, ve1] * corr[vs, ve2] # This is positive
                            localC_A_neg[vs] += tmp * np.absolute(pcorr)
                            denom_A_neg[vs] += tmp

    localC_A_pos /= denom_A_pos
    localC_A_neg /= denom_A_neg
    globalC_A_pos = np.average(localC_A_pos)
    globalC_A_neg = np.average(localC_A_neg)

    return globalC_A_pos, globalC_A_neg, localC_A_pos, localC_A_neg


# C^{cor,M}, i.e., clustering coef for correlation matrices based on partial mutual information
# Input
#     cov: covariance matrix, whose (i,j) element is stored in cov[nV*i+j]
# Output
#     Ci_M: vector of the local clustering coefficient
#     C_M: global clustering coefficient
def corrmat_C_M(cov):
    nV = len(cov) # number of nodes
    print(nV, 'nodes')
    corr = np.zeros((nV, nV))

    # calculate correlation matrix from covariance matrix
    for vs in range (nV):
        corr[vs, vs] = 1.0 # unused though
        for ve1 in range(vs):
            corr[vs, ve1] = cov[vs, ve1] / np.sqrt(cov[vs, vs] * cov[ve1, ve1])
            corr[ve1, vs] = corr[vs, ve1]

# The following is more concise but probably slower.
#
#    D_to_neg_half = np.diag(np.diag(cov) ** -0.5)
#    corr = np.dot(np.dot(D_to_neg_half, cov), D_to_neg_half)

    # initialization
    Ci_M = np.zeros(nV)
    denom = np.zeros(nV)

    # local clustering coef
    for vs in range(nV):
        for ve1 in range(nV):
            if ve1 != vs:
                for ve2 in range(ve1):
                    if ve2 != vs:
                        # partial mutual info between ve1 and ve2 controlling for vs
                        # X: ve1, Y:ve2, Z:vs. constant terms ignored. Gaussian signals assumed.
                        # Frenzel & Pompe, PRL (2007)
                        hXZ = 0.5 * np.log(1 - corr[ve1, vs]*corr[ve1, vs])
                        hYZ = 0.5 * np.log(1 - corr[ve2, vs]*corr[ve2, vs])
                        hZ = 0.0
                        hXYZ = 0.5 * np.log(1 - corr[ve2, vs]*corr[ve2, vs]
                            - corr[ve1, vs]*corr[ve1, vs] - corr[ve1, ve2]*corr[ve1, ve2]
                            + 2 * corr[ve1, ve2]*corr[ve1, vs]*corr[ve2, vs])
                        pmi = hXZ + hYZ - hZ - hXYZ
                        tmp = np.absolute(corr[vs, ve1] * corr[vs, ve2])
                        Ci_M[vs] += tmp * pmi
                        denom[vs] += tmp

#    Ci_M = Ci_M / (denom * 0.5 * (1 + np.log(2* math.pi)))
    Ci_M /= denom * 0.5 * (1 + np.log(2* math.pi))
    C_M = np.average(Ci_M)

    return C_M, Ci_M

# Clustering coef for correlation matrices, C^{cor,M}, but only using +++ triangles or --- triangles
#
# _pos: +++ triangles only
# _neg: --- triangles only
def C_corM_pmonly(cov):
    nV = len(cov) # number of nodes
    print(nV, 'nodes')
    corr = np.zeros((nV, nV))

    # calculate correlation matrix from covariance matrix
    for vs in range (nV):
        corr[vs, vs] = 1.0 # unused though
        for ve1 in range(vs):
            corr[vs, ve1] = cov[vs, ve1] / np.sqrt(cov[vs, vs] * cov[ve1, ve1])
            corr[ve1, vs] = corr[vs, ve1]

# The following is more concise but probably slower.
#
#    D_to_neg_half = np.diag(np.diag(cov) ** -0.5)
#    corr = np.dot(np.dot(D_to_neg_half, cov), D_to_neg_half)

    # initialization
    localC_pos = np.zeros(nV)
    localC_neg = np.zeros(nV)
    denom_pos = np.zeros(nV)
    denom_neg = np.zeros(nV)

    # local clustering coef
    for vs in range(nV):
        for ve1 in range(nV):
            if ve1 != vs:
                for ve2 in range(ve1):
                    if ve2 != vs:
                        # partial mutual info between ve1 and ve2 controlling for vs
                        # X: ve1, Y:ve2, Z:vs. constant terms ignored. Gaussian signals assumed.
                        # Frenzel & Pompe, PRL (2007)
                        hXZ = 0.5 * np.log(1 - corr[ve1, vs]*corr[ve1, vs])
                        hYZ = 0.5 * np.log(1 - corr[ve2, vs]*corr[ve2, vs])
                        hZ = 0.0
                        hXYZ = 0.5 * np.log(1 - corr[ve2, vs]*corr[ve2, vs]
                            - corr[ve1, vs]*corr[ve1, vs] - corr[ve1, ve2]*corr[ve1, ve2]
                            + 2 * corr[ve1, ve2]*corr[ve1, vs]*corr[ve2, vs])
                        pmi = hXZ + hYZ - hZ - hXYZ
                        if corr[vs, ve1] > 0 and corr[vs, ve2] > 0 and corr[ve1, ve2] > 0: # +++ triangle
                            tmp = corr[vs, ve1] * corr[vs, ve2] # This is positive.
                            localC_pos[vs] += tmp * pmi;
                            denom_pos[vs] += tmp;
                        elif corr[vs, ve1] < 0 and corr[vs, ve2] < 0 and corr[ve1, ve2] < 0: # --- triangle
                            tmp = corr[vs, ve1] * corr[vs, ve2] # This is positive.
                            localC_neg[vs] += tmp * pmi
                            denom_neg[vs] += tmp
                  
    localC_pos /= denom_pos * 0.5 * (1 + np.log(2*math.pi))
    localC_neg /= denom_neg * 0.5 * (1 + np.log(2*math.pi))
    globalC_pos = np.average(localC_pos)
    globalC_neg = np.average(localC_neg)

    return globalC_pos, globalC_neg, localC_pos, localC_neg
