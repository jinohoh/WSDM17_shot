#ifndef _MATRIX_
#define _MATRIX_

#include "Macros.hpp"
#include "Option.hpp"
#include "Util.hpp"

#include <omp.h>

#include <limits> // std::numeric_limits<PRECISION>
#include <cmath>

class DenseMatrix{
protected:
	uint szRow;
	uint szCol;
	DenseVector* mat;

public:
	PRECISION 	*ptrBulk;
	bool flagInit;
	DenseMatrix(): ptrBulk(nullptr), mat(NULL){
		szRow = 0;
		szCol = 0;
		flagInit = false;
	};
	DenseMatrix(uint row, uint col){
		init(row, col);
		flagInit = true;
	}
	~DenseMatrix(){}

	void cleanup(){
		szRow = 0;
		szCol = 0;
	}

	virtual DenseVector& getVec(uint id) = 0;
	void init(uint row, uint col){};
	uint getRowSize(){return szRow;}
	uint getColSize(){return szCol;}
};


class RowDenseMatrix : public DenseMatrix{
public:

	~RowDenseMatrix(){
		#pragma omp parallel for num_threads(Option::szThread)
		for (uint i = 0; i < szRow; i++)
			mat[i].ptr = nullptr;
		
		delete[] mat;
		free(ptrBulk);
	}

	void init(uint row, uint col, bool flagInitVal = true){
		cleanup();
		posix_memalign((void**)&ptrBulk, ALIGN_WIDTH, (long uint)row*col*sizeof(PRECISION));
		mat = new DenseVector[row];

		PRECISION* norm = new PRECISION[col];
		memset(ptrBulk, 0, (long uint)sizeof(PRECISION)*row*col);
		memset(norm, 0, (long long)sizeof(PRECISION)*col);

		#pragma omp parallel for num_threads(Option::szThread)
		for (uint i = 0; i < row; i++){
			mat[i].ptr = ptrBulk + (long long)i*col;
			if (flagInitVal){
				for (uint j = 0; j < col; j++){
					mat[i].ptr[j] = _drand();
					norm[j] += mat[i].ptr[j] * mat[i].ptr[j];
				}
			}
		}

		for (uint j = 0; j < col; j++){
			norm[j] = sqrt(norm[j]);
		}

		delete[] norm;

		szRow = row;
		szCol = col;
		flagInit = true;
	}

	void reinit(){
		memset(ptrBulk, 0, (long uint)sizeof(PRECISION)*szRow*szCol);
		return;
	}

	void norm(){
		double* proj = new double[szCol];
		memset(proj, 0, sizeof(double)*szCol);

		double denom;
		for (uint k = 0 ; k < szCol; k++){
			memset(proj, 0, sizeof(double)*szCol);

			#pragma omp parallel for num_threads(Option::szThread)
			for (uint j = k; j < szCol; j++)
				for (uint i = 0; i < szRow; i++){
					proj[j] += mat[i][k] * mat[i][j];
			}
			denom = sqrt(proj[k]);

			#pragma omp parallel for num_threads(Option::szThread)
			for (uint i = 0; i < szRow; i++){
				for (uint j = k+1; j < szCol; j++)
					if (0 != proj[k])
						mat[i][j] = mat[i][j] - proj[j]*mat[i][k]/proj[k];
				if (0 != denom)
					mat[i][k] = mat[i][k] / denom;
			}
		}
		
		delete[] proj;
	}

	void testNorm(){
		PRECISION *norms = new PRECISION[szCol];
		memset(norms, 0, sizeof(PRECISION)*szCol);
		for (uint i = 0; i < szRow; i++){
			for (uint j = 0; j < szCol; j++){
				norms[j] += mat[i][j] * mat[i][j];
			}
		}

		cout << "Normality check: ";
		for (uint j = 0; j < szCol; j++){
			norms[j] = sqrt(norms[j]);
			cout << norms[j] << " ";
		}

		cout << "\t" << "Orthonormal check: ";
		PRECISION proj = 0;
		for (uint p = 0; p < szCol; p++){
			for (uint j = p + 1; j < szCol; j++){
				proj = 0;
				for (uint i = 0; i < szRow; i++)
					proj += mat[i][p] * mat[i][j];
				
				cout << "(" << p << ", " << j << "): " << proj << "\t";
			}
		}
		cout << endl;
		delete[] norms;
	}

	DenseVector& getVec(uint rId){
		assert(rId < szRow);
		return mat[rId];
	}

	void updateColVec(uint cId, DenseVector& colvec){
		assert(cId < szCol);
		for (uint i = 0; i < szRow; i++)
			mat[i][cId] = colvec[i];
	}

	void genRandMat(int row, int col, float density = 1.0, bool flagSym = false, bool flagPD = false){
		if (flagSym){
			row = max(row, col);
			col = row;
		}
		init(row, col, false);

		PRECISION genVal = 0;
		for (uint i = 0; i < row; i++){
			for (uint j = 0; j < col; j++){
				if (_drand() >= density)
					continue;

				genVal = Option::initVar * _drand() - Option::initVar/2;
				mat[i][j] += genVal;

				if (flagSym)
					mat[j][i] += genVal;

				if (flagPD && (i == j))
					mat[i][i] += Option::epsilon;
			}
		}
	}
};

ostream& operator<<(ostream& os, RowDenseMatrix& rhs){

	typedef std::numeric_limits<PRECISION> lim;
	os.precision(32);

	uint szRow = rhs.getRowSize();
	uint szCol = rhs.getColSize();
	if (0 == szRow || 0 == szCol)
		return os;

	for (uint i = 0; i < szRow; i++){
		DenseVector& vec = rhs.getVec(i);
		os << scientific << vec[0];
		for (uint j = 1; j < szCol; j++)
			os << ", " << scientific << vec[j];
		os <<  endl;
	}

	return os;
}

#endif // #ifndef _MATRIX_
