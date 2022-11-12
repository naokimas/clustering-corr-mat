import numpy as np
import pandas as pd
import sys # to use input parameters and issue error

import weighted_clust_coef

df = pd.read_csv('cov2.txt', delimiter='\s+', header=None, index_col=None) # covariance matrix
nV = len(df) # number of rows
print(nV, "nodes")
cov = df.to_numpy() # convert to numpy array

# convert the covariance matrix into the edge list. If your data is already an edge list, you do not need this part
E = np.empty((0, 2), dtype=int) # edge list. Nodes E[i, 0] and E[i, 1] form an edge.
w = np.empty((0, 1)) # edge weight
for i in range(nV):
    for j in range(i):
        if cov[i, j] > 0:
            E = np.append(E, np.array([[i, j]]), axis=0)
            w = np.append(w, np.array([[cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])]]), axis=0) # edge weight
nE = len(E) # number of rows

C, localC = weighted_clust_coef.Cnet_zhang(nV, nE, E, w) # calculate the weighted clustering coefficient
print(C)
