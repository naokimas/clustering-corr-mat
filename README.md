# clustering-corr-mat
Clustering coefficients for correlation matrices. 
Please cite

Naoki Masuda, Michiko Sakaki, Takahiro Ezaki, Takamitsu Watanabe.<br>
Clustering coefficients for correlation networks.<br>
[Frontiers in Neuroinformatics, 12, 7 (2018).](https://doi.org/10.3389/fninf.2018.00007)

when using the code provided here.<br>
This paper proposes clustering coefficients for correlation matrix data.

## Code in Python

- **corrmat_clust_coef.py**: functions to calculate the local and global clustering coefficients for correlation matrices
- **test_corrmat_clust_coef.py**: an example script to use functions defined in corrmat_clust_coef.py
- **weighted_clust_coef.py**: functions to calculate the local and global weighted clustering coefficients for weighted networks. Other researchers proposed these methods in mid 2000s, and we examined them for correlation matrix data for comparison purposes in our paper.
- **test_weieghted_clust_coef.py**: an example script to use functions defined in weighted_clust_coef.py on a correlation matrix as input

Usage:

`python test_corrmat_clust_coef.py`

`python test_weighted_clust_coef.py`

To run test_corrmat_clust_coef.py or test_weighted_clust_coef.py, you need to prepare a correlation matrix data file as a matrix and place it in the same folder. Change the input file name in test_corrmat_clust_coef.py or test_weighted_clust_coef.py accordingly.

## Code in C++

**corrmat_cc_group.cp** provides functions to calculate local and global clustering coefficients for correlation matrices and some auxiliary functions.

Function corrmat_C_PHA calculates C^{cor,P}, C^{cor,H}, C^{cor,A} and their local variants (i.e. C_i's)

Function corrmat_C_M calculates C^{cor,M} and its local variant

Example usage is shown after the definition of the two functions.
