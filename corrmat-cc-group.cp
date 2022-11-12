/* Local and global clustering coefficients for correlation matrices.
   Edge list must not be duplicated.
    Output:
        C-P, C-H, C-A, C-M, s, s-P:
            i-th row: i-th ROI
            j-th column: j-th participant
        si-vs-ci:
            1st column: participant ID (starting from 1)
            2nd column: s_i
            3rd column: s^+_i
            4th column: C^{cor,P}_i
            5th column: C^{cor,H}_i
            6th column: C^{cor,A}_i
            7th column: C^{cor,M}_i
*/

#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cmath> // sqrt

#include "../../social/lib/stat/mt19937ar.c"
#include "tools-cc-corrmat.cp"

int main (int argc, char **argv) {

	if (argc != 2 && argc != 3) {
		cerr << "Usage: corrnet-cc-group.out data mode(optional)" << endl;
		exit(8);
	}

    int nP = 28; // # participants
    int nV = 30; // # nodes
    int data = atoi(argv[1]);
    std::string basename;
    if (data == 0) { // younger adults
        basename = "YO/Y";
        cerr << "younger adults (30 ROIs)" << endl;
    } else if (data == 1) { // older adults
        basename = "YO/O";
        cerr << "older adults (30 ROIs)" << endl;
    } else if (data == 2) { // younger adults
        basename = "YO264/Y";
        nV = 264;
        cerr << "younger adults (264 ROIs)" << endl;
    } else if (data == 3) { // older adults
        basename = "YO264/O";
        nV = 264;
        cerr << "older adults (264 ROIs)" << endl;
    } else if (data == 4) { // typically developing adults (as opposed to ASD)
        basename = "TDASD/TD";
        cerr << "TD" << endl;
        nP = 26;
    } else if (data == 5) { // ASD individuals
        basename = "TDASD/ASD";
        cerr << "ASD" << endl;
        nP = 26;
    } else if (data == 6) { // continuous age
        basename = "NKI138/cov";
        cerr << "younger and older (continuous age)" << endl;
        nP = 138;
    } else if (data == 8) { // continuous age. motor + auditory + visual networks from Power
        basename = "NKI138/SM";
        cerr << "continuous age, motor + auditory + visual" << endl;
        nP = 138;
        nV = 78;
    } else if (data == 9) { // continuous age, global signal = gray matter only 
        basename = "NKI138/covGM";
        cerr << "younger and older (continuous age), global signal = GM only" << endl;
        nP = 138;
    } else {
        cerr << "data must be 0 (younger, 30 ROIs), 1 (older, 30 ROIs), 2 (younger, 264 ROIs), 3 (older, 264 ROIs), 4 (TD, 30 ROIs), 5 (ASD, 30 ROIs), or 6 (younger and older, continuous age, 30 ROIs), 8 (continuous age, motor + auditory + visual networks from Power, 78 ROIs), 9 (younger and older, continuous age, 30 ROIs, global signal = GM only)" << endl;
        exit(8);
    }

    int mode;
    if (argc == 2)
        mode = 0; // use original data
    else { // argc == 3    
        mode = atoi(argv[2]);
        if (mode >= 3 || mode <= -1) {
            cerr << "mode must be 0 (use original data), 1 (HQS), or 2 (white noise)" << endl;
            exit(8);
        }
    }

    int i, j, ii, ip; // counters
    double crosscorr; // working var
    double cov[nV*nV]; // covariance matrix

    double cov_out[nV*nV]; // working var; used only when mode == 1 (i.e. HQS) 
    if (mode==1)
        init_genrand(time(NULL));

    double Ci_P[nV*nP], Ci_H[nV*nP], Ci_A[nV*nP], Ci_M[nV*nP], Ci_P_tmp[nV], Ci_H_tmp[nV], Ci_A_tmp[nV], Ci_M_tmp[nV]; // local clustering coef
    double globalC_P[nP], globalC_H[nP], globalC_A[nP], globalC_M[nP], globalC_P_tmp, globalC_H_tmp, globalC_A_tmp, globalC_M_tmp; // global clustering coef
	double s[nV*nP], s_P[nV*nP], global_s[nP], global_s_P[nP]; // (signed) node strength

    int successful_HQS = 0; // becomes -1 if unsuccessful due to non-positive avg off-diagonal
    ofstream fout("si-vs-ci");

    for (ip=0 ; ip<nP ; ip++) { // ip-th participant

        // read the covaraince matrix from the input file
        read_to_cov(basename + std::to_string(ip+1) + ".txt", nV, cov);

        if (mode == 1) { // HQS
            successful_HQS = H_Q_S(cov, cov_out, nV); // = 0 if successful, and = -1 if unsuccessful
            if (successful_HQS == -1)
                cerr << "HQS algorithm failed for participant " << ip+1 << " due to non-positive avg off-diagonal" << endl;
            for (i=0 ; i<nV ; i++)
                for (j=0 ; j<nV ; j++)
                    cov[i*nV+j] = cov_out[i*nV+j];
        } else if (mode == 2) { // white noise
            white_cov(cov, nV, 0.0, 200);
        }

        // clustering coef for correlation matrix
        corrmat_C_PHA(cov, nV, Ci_P_tmp, Ci_H_tmp, Ci_A_tmp, &globalC_P_tmp, &globalC_H_tmp, &globalC_A_tmp);

        corrmat_C_M(cov, nV, Ci_M_tmp, &globalC_M_tmp);

        for (i=0 ; i<nV ; i++) {
            Ci_P[ip*nV+i] = Ci_P_tmp[i];
            Ci_H[ip*nV+i] = Ci_H_tmp[i];
            Ci_A[ip*nV+i] = Ci_A_tmp[i];
            Ci_M[ip*nV+i] = Ci_M_tmp[i];
        }
        globalC_P[ip] = globalC_P_tmp;
        globalC_H[ip] = globalC_H_tmp;
        globalC_A[ip] = globalC_A_tmp;
        globalC_M[ip] = globalC_M_tmp;

        // node strength
        for (i=0 ; i<nV ; i++)
            s[ip*nV+i] = s_P[ip*nV+i] = 0.0;
        for (i=0 ; i<nV ; i++)
            for (j=i+1 ; j<nV ; j++) {
                crosscorr = cov[i*nV+j]/sqrt(cov[i*nV+i])/sqrt(cov[j*nV+j]); // cross-correlation
                s[ip*nV+i] += crosscorr;
                s[ip*nV+j] += crosscorr;
                if (crosscorr>0) {
                    s_P[ip*nV+i] += crosscorr;
                    s_P[ip*nV+j] += crosscorr;
                }
            }
        global_s[ip] = global_s_P[ip] = 0.0;
        for (i=0 ; i<nV ; i++) {
            s[ip*nV+i] /= nV-1;
            s_P[ip*nV+i] /= nV-1;
            global_s[ip] += s[ip*nV+i];
            global_s_P[ip] += s_P[ip*nV+i];
        }
        global_s[ip] /= nV;
        global_s_P[ip] /= nV;

        if (successful_HQS == 0) { // successful
            for (i=0 ; i<nV ; i++) {
               fout << ((data != 6)? ip+1 + nP*data : ip+1) << " " << s[ip*nV + i] << " " << s_P[ip*nV+i] << " " << Ci_P[ip*nV + i] << " " << Ci_H[ip*nV + i] << " " << Ci_A[ip*nV + i] << " " << Ci_M[ip*nV+i] << endl;
            }
        } else { // unsuccessful
            global_s_P[ip] = NAN; // force nan
        }
    } // all participants done

    // output local quantities in the format useful for statistical analysis
    ofstream fout_C_P("C-P");
    ofstream fout_C_H("C-H");
    ofstream fout_C_A("C-A");
    ofstream fout_C_M("C-M");
    ofstream fout_s("s");
    ofstream fout_s_P("s-P");
    for (i=0 ; i<nV ; i++) {
        for (ip=0 ; ip<nP ; ip++) {
            if (ip>0) {
                fout_C_P << " ";
                fout_C_H << " ";
                fout_C_A << " ";
                fout_C_M << " ";
                fout_s << " ";
                fout_s_P << " ";
            }
            fout_C_P << Ci_P[ip*nV + i];
            fout_C_H << Ci_H[ip*nV + i];
            fout_C_A << Ci_A[ip*nV + i];
            fout_C_M << Ci_M[ip*nV + i];
            fout_s << s[ip*nV + i];
            fout_s_P << s_P[ip*nV + i];
        }
        fout_C_P << endl;
        fout_C_H << endl;
        fout_C_A << endl;
        fout_C_M << endl;
        fout_s << endl;
        fout_s_P << endl;
    }
    fout_C_P.close();
    fout_C_H.close();
    fout_C_A.close();
    fout_C_M.close();
    fout_s.close();
    fout_s_P.close();

    // output global quantities
    ofstream fout_global("corrmat-c");
    for (ip=0 ; ip<nP ; ip++) {
        fout_global << ((data !=  6 && data != 9)? ip+1 + nP*data : ip+1) << " " << globalC_P[ip] << " " << globalC_H[ip] << " " << globalC_A[ip] << " " << globalC_M[ip] << " " << global_s[ip] << " " << global_s_P[ip] << endl;
    }
    fout_global.close();

    fout.close();
    return 0;
}