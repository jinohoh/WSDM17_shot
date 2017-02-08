#ifndef _REVEIG_
#define _REVEIG_

#include "Macros.hpp"
#include <iostream>	// cout <<
#include <fstream> 	// fstream open()
#include <cstdio> 	// sprintf()
#include "arrssym.h" // ARrcSymStdEig: ARpack code for reverse communication symmetric standard eigen problem.

using namespace std;
typedef ARrcSymStdEig<PRECISION> ArSolver;

inline void testNaiveRevEig(){

	int eigenSize = 2;
	int colSize = 10;

	PRECISION **mat = new PRECISION*[colSize];
	mat[0] = new PRECISION[colSize*colSize];
	for (uint i = 1; i < colSize; i++)
		mat[i] = mat[0] + i * colSize;

	for (uint i = 0; i < colSize; i++)
		for (uint j = 0; j < colSize; j++)
			mat[i][j] = (i + j)*colSize + (i + j);

	ARrcSymStdEig<PRECISION>  prob(colSize, eigenSize, "LM");
	while( !prob.ArnoldiBasisFound()){
		if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) { 
			PRECISION* wold = prob.GetVector();
			PRECISION* wnew = prob.PutVector();
			// Multiplicate
			for (uint i = 0; i < colSize; i++){
				wnew[i] = 0;
				for (uint j = 0; j < colSize; j++)
					wnew[i] += mat[i][j]*wold[j];
			}
		}
		prob.TakeStep();
	}

	// Finalize Eigenvectors 
	prob.FindEigenvectors();

	// Codes from Solution().
	uint  i, n, nconv, mode;
	n     = prob.GetN();
	nconv = prob.ConvergedEigenvalues();
	mode  = prob.GetMode();
	
	cout << endl << endl << "Testing ARPACK++ class ARrcSymStdEig" << endl;
	cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << endl;
	switch (mode) {
	case 1:
	  cout << "Regular mode" << endl << endl;
	  break;
	case 3: 
	  cout << "Shift and invert mode" << endl << endl;
	}
	
	cout << "Dimension of the system            : " << n              << endl;
	cout << "Number of 'requested' eigenvalues  : " << prob.GetNev()  << endl;
	cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
	cout << "Number of Arnoldi vectors generated: " << prob.GetNcv()  << endl;
	cout << "Number of iterations taken         : " << prob.GetIter() << endl;
	cout << endl;
	
	if (prob.EigenvaluesFound()) {
	  cout << "Eigenvalues:" << endl;
	  for (i=0; i<nconv; i++) {
	    cout << "  lambda[" << (i+1) << "]: " << prob.Eigenvalue(i) << endl;
	  }
	  cout << endl;
	}

	// update Eigenvector 
	cout << "Converged Eigenvector" << endl;
	for (uint i = 0; i < eigenSize ; i++){
		cout << i << "th Eigenvector\t";
		for (uint j = 0; j < colSize; j++)
			cout << prob.Eigenvector(i, j) << "\t";
		cout << endl;
	}
}

inline void testRandom(const char* filename){

	Option::initVar = 2.;
	Option::epsilon = 1.;
	char temp[256];
	char temp2[16];
	int matNum = 0;

	RowDenseMatrix randMat;
	fstream outFP;

	for (int matNum = 0; matNum < 1; matNum++){

		// Regular matrix generation
		randMat.genRandMat(30, 5);
		strcpy(temp, filename);
		sprintf(temp2, ".reg.%d", matNum);
		strcat(temp, temp2);
		cout << "New filename: " <<  temp << endl;
		outFP.open(temp, fstream::out | fstream::trunc);
		outFP << randMat;
		outFP.close();
		cout << randMat;

		// psd matrix generation
		randMat.genRandMat(30, 5, 1, true, true);
		strcpy(temp, filename);
		sprintf(temp2, ".psd.%d", matNum);
		strcat(temp, temp2);
		cout << "New filename: " <<  temp << endl;
		outFP.open(temp, fstream::out | fstream::trunc);
		outFP << randMat;
		outFP.close();
		cout << randMat;
	}
}



#endif // #ifndef _REVEIG_
