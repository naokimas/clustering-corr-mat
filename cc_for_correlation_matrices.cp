/* Local and global clustering coefficients for correlation matrices.

    Function corrmat_C_PHA calculates C^{cor,P}, C^{cor,H}, C^{cor,A} and their local variants (i.e. C_i's)

    Function corrmat_C_M calculates C^{cor,M} and its local variant

    Example usage is shown after the two functions. */


/* clustering coef for correlation matrices based on the partial correlation coefficient
Input
    cov: covariance matrix, whose (i,j) element is stored in cov[nV*i+j]
    nV: number of nodesa
Output
    Ci_{type}: local clustering coefficient. An array of size nV should be supplied for storing the results.
    C_{type}: global clustering coefficient. The pointer to a variable should be supplied for storing the results.
    type = P, H, or A
    P: C^{cor,P}. "P" for positive. Consider only positive edges.
    H: C^{cor,H}. "H" for Heider. Consider all triangles and use the Heider's balance theory criterion.
    A: C^{cor,A}. "A" for absolute value. Use the absolute value of the partial correlation coefficient.
*/
int corrmat_C_PHA(double *cov, int nV, double *Ci_P, double *Ci_H, double *Ci_A, double *C_P, double *C_H, double *C_A) {

    int vs,ve1,ve2; // counters
    double corr[nV*nV], tmp, pcorr; // working vars
    double denom_P[nV], denom_HA[nV]; // denominators

    // calculate correlation matrix from covariance matrix
    for (vs=0 ; vs<nV ; vs++) {
        corr[vs*nV + vs] = 1.0; // unused though
        for (ve1=0 ; ve1<vs ; ve1++) {
            corr[vs*nV + ve1] = cov[vs*nV + ve1] / sqrt(cov[vs*nV + vs] * cov[ve1*nV + ve1]);
            corr[ve1*nV + vs] = corr[vs*nV + ve1];
        }
    } 

    // initialization
    for (vs=0 ; vs<nV ; vs++)
        Ci_P[vs] = Ci_H[vs] = Ci_A[vs] = denom_P[vs] = denom_HA[vs] = 0.0;

    // local clustering coef
    for (vs=0 ; vs<nV ; vs++)
        for (ve1=0 ; ve1<nV ; ve1++)
            if (ve1 != vs) {
                for (ve2=0 ; ve2<ve1 ; ve2++) {
                    if (ve2 != vs) {
                        // partial corr between ve1 and ve2 controlling for vs
                        pcorr = (corr[ve1*nV+ve2] - corr[vs*nV+ve1] * corr[vs*nV+ve2]) /
                            sqrt( (1 - corr[vs*nV+ve1] * corr[vs*nV+ve1]) *
                            (1 - corr[vs*nV+ve2] * corr[vs*nV+ve2]) );
                        tmp = corr[vs*nV+ve1] * corr[vs*nV+ve2];
                        Ci_H[vs] += tmp * pcorr; // Heider
                        Ci_A[vs] += (tmp * pcorr > 0) ?  tmp * pcorr : - tmp * pcorr; // absolute value
                        denom_HA[vs] += (tmp>0) ? tmp : -tmp; // add abs(tmp)
                        if (corr[vs*nV+ve1] > 0 && corr[vs*nV+ve2] > 0) { // (vs,ve1) and (vs,ve2) are both positive edges
                            Ci_P[vs] += tmp * pcorr;
                            denom_P[vs] += tmp;
                        }
                    }
                }
            }

    *C_P = *C_H = *C_A = 0.0; // initialization
    for (vs=0 ; vs<nV ; vs++) {
        Ci_P[vs] /= denom_P[vs];
        Ci_H[vs] /= denom_HA[vs];
        Ci_A[vs] /= denom_HA[vs];
        *C_P += Ci_P[vs];
        *C_H += Ci_H[vs];
        *C_A += Ci_A[vs];
    }
    *C_P /= nV;
    *C_H /= nV;
    *C_A /= nV;
    return 0;
}


/* C^{cor,M}, i.e., clustering coef for correlation matrices based on partial mutual information
Input
    cov: covariance matrix, whose (i,j) element is stored in cov[nV*i+j]
    nV: number of nodesa
Output
    Ci_M: local clustering coefficient. An array of size nV should be supplied for storing the results.
    C_M: global clustering coefficient. The pointer to a variable should be supplied for storing the results.
*/
int corrmat_C_M(double *cov, int nV, double *Ci_M, double *C_M) {
    int vs,ve1,ve2; // counters
    double corr[nV*nV], tmp, pcorr; // working vars
    double denom[nV]; // denominators

    // calculate correlation matrix from covariance matrix
    for (vs=0 ; vs<nV ; vs++) {
        corr[vs*nV + vs] = 1.0; // unused though
        for (ve1=0 ; ve1<vs ; ve1++) {
            corr[vs*nV + ve1] = cov[vs*nV + ve1] / sqrt(cov[vs*nV + vs] * cov[ve1*nV + ve1]);
            corr[ve1*nV + vs] = corr[vs*nV + ve1];
        }
    } 

    // initialization
    for (vs=0 ; vs<nV ; vs++)
        Ci_M[vs] = denom[vs] = 0.0;
    double hZ, hXZ, hYZ, hXYZ, pmi; // working vars

    // local clustering coef
    for (vs=0 ; vs<nV ; vs++)
        for (ve1=0 ; ve1<nV ; ve1++)
            if (ve1 != vs) {
                for (ve2=0 ; ve2<ve1 ; ve2++) {
                    if (ve2 != vs) {
                        // partial mutual info between ve1 and ve2 controlling for vs
                        // X: ve1, Y:ve2, Z:vs. constant terms ignored. Gaussian signals assumed.
                        // Frenzel & Pompe, PRL (2007)
                        hXZ = 0.5 * log(1 - corr[ve1*nV+vs]*corr[ve1*nV+vs]);
                        hYZ = 0.5 * log(1 - corr[ve2*nV+vs]*corr[ve2*nV+vs]);
                        hZ = 0.0;
                        hXYZ = 0.5 * log(1 - corr[ve2*nV+vs]*corr[ve2*nV+vs]
                            - corr[ve1*nV+vs]*corr[ve1*nV+vs] - corr[ve1*nV+ve2]*corr[ve1*nV+ve2]
                            + 2 * corr[ve1*nV+ve2]*corr[ve1*nV+vs]*corr[ve2*nV+vs]);
                        pmi = hXZ + hYZ - hZ - hXYZ;
                        tmp = corr[vs*nV+ve1] * corr[vs*nV+ve2];
                        if (tmp<0)
                            tmp = -tmp; // now tmp >= 0
                        Ci_M[vs] += tmp * pmi;
                        denom[vs] += tmp;
                    }
                }
            }

    *C_M = 0.0;
    for (vs=0 ; vs<nV ; vs++) {
        Ci_M[vs] /= denom[vs] * 0.5 * (1 + log(2*M_PI));
        *C_M += Ci_M[vs];
    }
    *C_M /= nV;
    return 0;
}

/* Example usage */
#ifdef MAIN
int main (int argc, char **argv) {

    int nV = 100; // # nodes
    double cov[nV*nV]; // covariance matrix

    /* set cov[] */

    double Ci_P[nV], Ci_H[nV], Ci_A[nV], Ci_M[nV]; // local clustering coefficients
    double globalC_P, globalC_H, globalC_A, globalC_M; // global clustering coefficients

    // Calculate C^{cor,A}, C^{cor,P}, C^{cor,H} and their local variants
    corrmat_C_PHA(cov, nV, Ci_P, Ci_H, Ci_A, &globalC_P, &globalC_H, &globalC_A);

    // Calculate C^{cor,M} and its local variant
    corrmat_C_M(cov, nV, Ci_M, &globalC_M);

    /* output or further processing */

    return 0;
}
#endif