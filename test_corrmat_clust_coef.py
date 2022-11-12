import numpy as np
import pandas as pd
import sys # to use input parameters and issue error

import corrmat_clust_coef

df = pd.read_csv('cov1.txt', delimiter='\s+', header=None, index_col=None) # read covariance matrix data
nV = len(df) # number of rows
print(nV, "nodes")
cov = df.to_numpy() # convert to numpy array. cov is an nV \times nV matrix
C_M, Ci_M = corrmat_clust_coef.corrmat_C_M(cov) # calculate the clustering coefficient, the one based on the partial mutual information, for the correlation matrix
print(C_M, Ci_M)
