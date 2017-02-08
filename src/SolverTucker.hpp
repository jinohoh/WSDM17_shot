#ifndef _ENGINE_TUCKER_FACT_
#define _ENGINE_TUCKER_FACT_

#include <iostream>
#include <cmath>
#include <omp.h>
#include "TensorFile.hpp"
#include "Matrix.hpp"
#include "RevEigSolver.hpp"
#include "Option.hpp"

using namespace std;

template <typename Tuple>
uint TensorIndexWithMultipleFile<Tuple>::Scanner::sCurMode;
template <typename Tuple>
char TensorIndexWithMultipleFile<Tuple>::Scanner::snScanner;

template <typename Tuple>
class SolverTuckerLS 
{
protected:
	uint mode;
	uint nconv;
	uint lenV;
	ArSolver 							*arProb; 

	SecureDenseVector					vOldAVX;
	SecureDenseVector					vNewAVX;
	PRECISION							*eigValues;

	// Variables created by ohter classes.
	TensorIndexWithMultipleFile<Tuple>	*ptrIndex;
	RowDenseMatrix						*factMat;

	uint smallit;
public:
	SolverTuckerLS(uint inMode, TensorIndexWithMultipleFile<Tuple>* inPtrIndex, RowDenseMatrix *inFactMat)
					: mode(inMode), arProb(nullptr), eigValues(nullptr), ptrIndex(inPtrIndex), factMat(inFactMat) {}

	virtual ~SolverTuckerLS()
	{
		delete arProb;
		delete[] eigValues;
	}

	void init()
	{
		lenV = this->getLenV(); 
		arProb = new ArSolver(lenV, Option::szRank[mode], "LM"); 
		arProb->TakeStep();
		vOldAVX.resize(lenV);
		vNewAVX.resize(lenV);
		eigValues = new PRECISION[Option::szRank[mode]];
		memset(eigValues, 0, Option::szRank[mode]*sizeof(PRECISION));
	}

	double doIter()
	{
		lenV = this->getLenV(); 
		smallit = 0;

		arProb->Restart();
		arProb->TakeStep();
		while (!arProb->ArnoldiBasisFound())
		{
			PRECISION* vOld = arProb->GetVector();
			PRECISION* vNew = arProb->PutVector();
			memcpy(vOldAVX.ptr, vOld, sizeof(PRECISION)*lenV);
			this->updatev(vOldAVX.ptr, vNewAVX.ptr);
			memcpy(vNew, vNewAVX.ptr, sizeof(PRECISION)*lenV);
			arProb->TakeStep();
			smallit++;
		}

		nconv = arProb->FindEigenvectors();

		// Compute the amount of change in Eigenvalues.
		double squareLossSum = 0;
		double squareSum = 0;
		double diff;
		for (uint i = 0 ; i < Option::szRank[mode]; i++)
		{
			diff = eigValues[i] - arProb->Eigenvalue(i);
			eigValues[i] = arProb->Eigenvalue(i);
			squareLossSum += diff*diff;
			squareSum += (double)eigValues[i] * eigValues[i];
		}
		if (squareSum != 0)
			squareLossSum = sqrt(squareLossSum)/sqrt(squareSum);
		else
			squareLossSum = sqrt(squareLossSum);
		// Store solution into factMat
		this->storeSol();

		return squareLossSum;
	}

	virtual uint getLenV()
	{
		return 1;
	}
	virtual void updatev(PRECISION *vOld, PRECISION *vNew){}
	virtual void storeSol(){};
};

template <typename Tuple>
class SolverTransLS : public SolverTuckerLS<Tuple>
{
	SecureDenseVector 					*segVNew;
	SecureDenseVector					*tempMat;
	bool								initFlag;

	typedef typename TensorIndexWithMultipleFile<Tuple>::Scanner Scanner;
	using SolverTuckerLS<Tuple>::mode;
	using SolverTuckerLS<Tuple>::nconv;
	using SolverTuckerLS<Tuple>::lenV;
	using SolverTuckerLS<Tuple>::arProb;
	using SolverTuckerLS<Tuple>::ptrIndex;
	using SolverTuckerLS<Tuple>::factMat;
	using SolverTuckerLS<Tuple>::smallit;
public:

	SolverTransLS(uint inMode, TensorIndexWithMultipleFile<Tuple>* inPtrIndex, RowDenseMatrix *inFactMat)
					: SolverTuckerLS<Tuple>(inMode, inPtrIndex, inFactMat), initFlag(true)
	{
		tempMat = new SecureDenseVector[Option::szRank[mode]];
		segVNew = new SecureDenseVector[Option::szThread];
	}

	~SolverTransLS()
	{
		delete[] tempMat;
		delete[] segVNew;
	}

	uint getLenV()
	{
		this->lenV = 1;
		for (uint i = 0; i < Option::szMode; i++)
			if (i != mode)
				lenV = lenV * Option::szRank[i];
		return lenV;
	}

	void updatev(PRECISION *vOld, PRECISION *vNew)
	{
		// Update an arnoldi vector.
		Scanner::setMode(mode);
		uint szRow = 0;
		uint nnz = 0;
		#pragma omp parallel for num_threads(Option::szThread)
		for (int tid = 0; tid < Option::szThread; tid++)
		{
			segVNew[tid].resize(lenV, 0);
			Tuple temp;
			int lastRowID = -1;
			PRECISION row_dot_wold;
			SecureDenseVector rowVec(lenV);
			SecureDenseVector tempVec(lenV);
			Scanner scan(ptrIndex);
			scan.rewind();
			while(scan.getNext(temp))
			{
				nnz++;
				if (lastRowID != temp.dim[mode])
				{
					if (-1 != lastRowID)
					{	// When a rowVec is finalized.
						row_dot_wold = innerProd(rowVec, vOld);
						segVNew[tid].acc(row_dot_wold, rowVec);
						szRow++;
					}

					// Init variables for next it.
					lastRowID = temp.dim[mode];
					rowVec.resize(lenV, 0);
				}

				tempVec.resize(1, temp.value);
				for (uint mIter = 0; mIter < Option::szMode; mIter++)
					if (mIter != mode)
						tempVec.outerProd(factMat[mIter].getVec(temp.dim[mIter]), Option::szRank[mIter]);
				rowVec.acc(tempVec);
			}
			if (-1 != lastRowID){
				row_dot_wold = innerProd(rowVec, vOld);
				segVNew[tid].acc(row_dot_wold, rowVec);
				szRow++;
			}
		}

		// Aggregating results from multiple threads.
		memset(vNew, 0, sizeof(PRECISION)*lenV);
		for (int tid = 0; tid < Option::szThread; tid++)
			for (int j = 0; j < lenV; j++)
				vNew[j] += segVNew[tid][j];


		return;
	}

	void storeSol(){

		uint maxEigIdx = nconv - 1;
		uint szEigVec = min(nconv, Option::szRank[mode]);

		//cout << "Mode: " << mode << "\t(" << lenV << ", " << maxEigIdx+1 << "), " << arProb->GetIter() << ", " << smallit << " steps";
		//for (uint j = 0; j < szEigVec; j++){
		//	cout << "\t" << arProb->Eigenvalue(maxEigIdx - j);
		//	double value = arProb->Eigenvalue(maxEigIdx - j);
		//	if (value < 1E-35){
		//		szEigVec = j;
		//	}
		//}
		//cout << endl;

		if (initFlag){
			factMat[mode].reinit();
			initFlag = false;
		}

		for (uint j = 0; j < szEigVec; j++){
			tempMat[j].resize(lenV, 0);
			for (uint i = 0; i < lenV; i++)
				tempMat[j][i] = arProb->Eigenvector(maxEigIdx - j, i);
		}


		/*
		 * computing eigen of original prob. based on that of trans. requires one extra scanning.
		 */
		Scanner::setMode(mode);
		#pragma omp parallel for num_threads(Option::szThread)
		for (uint tid = 0; tid < Option::szThread; tid++){
			Tuple temp;
			int lastRowID = -1;
			SecureDenseVector rowVec(lenV);
			SecureDenseVector tempVec(lenV);
			Scanner scan(ptrIndex);
			scan.rewind();
			while(scan.getNext(temp)){
				if (-1 == lastRowID)
					lastRowID = temp.dim[mode];
				else if (lastRowID != temp.dim[mode]){
					// When a rowVec is finalized.
					for (uint j = 0; j < szEigVec; j++){
						factMat[mode].getVec(lastRowID)[j] = innerProd(rowVec, tempMat[j]);
					}
					lastRowID = temp.dim[mode];
					rowVec.resize(lenV, 0);
				}

				tempVec.resize(1, temp.value);
				for (uint mIter = 0; mIter < Option::szMode; mIter++)
					if (mIter != mode)
						tempVec.outerProd(factMat[mIter].getVec(temp.dim[mIter]), Option::szRank[mIter]);
				rowVec.acc(tempVec);
			}
			if (-1 != lastRowID){
				for (uint j = 0; j < szEigVec; j++){
					factMat[mode].getVec(lastRowID)[j] = innerProd(rowVec, tempMat[j]);
				}
			}
		}

		if (Option::flagFlip){
			#pragma omp parallel for num_threads(Option::szThread)
			for (uint j = 0; j < szEigVec; j++){
				PRECISION maxMagnitude = 0;
				for (uint i = 0; i < Option::szDim[mode]; i++){
					if (fabs(factMat[mode].getVec(i)[j]) > fabs(maxMagnitude))
						maxMagnitude = factMat[mode].getVec(i)[j];
				}
				if (maxMagnitude < 0)
					for (uint i = 0; i < Option::szDim[mode]; i++)
						factMat[mode].getVec(i)[j] = -factMat[mode].getVec(i)[j];
			}
		}

		factMat[mode].norm();
	}

};

template <typename Tuple>
class SolverRowLS: public SolverTuckerLS<Tuple>
{

	SecureDenseVector *segSumVec;

	typedef typename TensorIndexWithMultipleFile<Tuple>::Scanner Scanner;
	using SolverTuckerLS<Tuple>::mode;
	using SolverTuckerLS<Tuple>::nconv;
	using SolverTuckerLS<Tuple>::lenV;
	using SolverTuckerLS<Tuple>::arProb;
	using SolverTuckerLS<Tuple>::ptrIndex;
	using SolverTuckerLS<Tuple>::factMat;
	using SolverTuckerLS<Tuple>::smallit;

public:
	SolverRowLS(uint inMode, TensorIndexWithMultipleFile<Tuple>* inPtrIndex, RowDenseMatrix *inFactMat)
					: SolverTuckerLS<Tuple>(inMode, inPtrIndex, inFactMat) 
	{
		segSumVec = new SecureDenseVector[Option::szThread];
	}

	~SolverRowLS(){
		delete[] segSumVec;
	}

	uint getLenV(){
		return Option::szDim[mode];
	}

	void updatev(PRECISION *vOld, PRECISION *vNew){
		uint lenRowVec = 1;
		for (uint i = 0; i < Option::szMode; i++)
			if (i != mode)
				lenRowVec = lenRowVec * Option::szRank[i];

		/*
		 * 1. Computing sumVec = Y'v
		 */
		Scanner::setMode(mode);

		SecureDenseVector sumVec;
		uint nnz = 0;
		#pragma omp parallel for num_threads(Option::szThread)
		for (int tid = 0; tid < Option::szThread; tid++){
			segSumVec[tid].resize(lenRowVec, 0);
			Tuple temp;

			SecureDenseVector tempVec(lenRowVec);
			Scanner scan(ptrIndex);
			scan.rewind();
			while(scan.getNext(temp)){
				// Computing a fraction of row vector (outer product of factor vectors)
				tempVec.resize(1, temp.value*vOld[temp.dim[mode]]);
				for (int mIter = 0 ; mIter < Option::szMode; mIter++)
					if (mIter != mode)
						tempVec.outerProd(factMat[mIter].getVec(temp.dim[mIter]), Option::szRank[mIter]);

				segSumVec[tid].acc(tempVec);
				nnz++;
			}
		}
		sumVec.resize(lenRowVec, 0);
		for (int tid = 0; tid < Option::szThread; tid++)
			sumVec.acc(segSumVec[tid]);

		/*
		 * 2. Computing vNew = Y * sumVec
		 */
		Scanner::setMode(mode);
		memset(vNew, 0, sizeof(PRECISION)*lenV);
		#pragma omp parallel for num_threads(Option::szThread)
		for (int tid = 0; tid < Option::szThread; tid++){
			Tuple temp;
			SecureDenseVector tempVec(lenRowVec);
			Scanner scan(ptrIndex);
			scan.rewind();

			// Updating value while scanning a data.
			while(scan.getNext(temp)){
				tempVec.resize(1, temp.value);
				for (uint mIter = 0; mIter < Option::szMode; mIter++)
					if (mIter != mode)
						tempVec.outerProd(factMat[mIter].getVec(temp.dim[mIter]), Option::szRank[mIter]);
				vNew[temp.dim[mode]] += innerProd(tempVec, sumVec);
			}
		}

		return;
	}

	void storeSol(){
		uint maxEigIdx = min(nconv-1, Option::szRank[mode]-1);
		uint szEigVec = min(nconv, Option::szRank[mode]);

		PRECISION maxMagnitude;
		for (uint j = 0; j < szEigVec; j++){
			maxMagnitude = 0;
			#pragma omp parallel for num_threads(Option::szThread)
			for (uint i = 0; i < lenV; i++){
				factMat[mode].getVec(i)[j] = arProb->Eigenvector(maxEigIdx - j, i);
				if (fabs(factMat[mode].getVec(i)[j]) > fabs(maxMagnitude))
					maxMagnitude = factMat[mode].getVec(i)[j];
			}

			if (Option::flagFlip && maxMagnitude < 0)
				#pragma omp parallel for num_threads(Option::szThread)
				for (uint i = 0; i < lenV; i++)
					factMat[mode].getVec(i)[j] = -factMat[mode].getVec(i)[j];
		}

		/*
		 * Printing results
		 */
		//cout << "Mode: " << mode << "\t(" << lenV << ", " << maxEigIdx+1 << "), " << arProb->GetIter() << ", " << smallit << " steps";
		//for (uint j = 0; j < szEigVec; j++){
		//	cout << "\t" << arProb->Eigenvalue(maxEigIdx - j);
		//}
		//cout << endl;

		return;
	}

};


template<typename Tuple>
class SolverTuckerFact
{
protected:
	// instance FactorModel
	typedef typename TensorIndexWithMultipleFile<Tuple>::Scanner Scanner;
	TensorIndexWithMultipleFile<Tuple>	*ptrIndex;
	RowDenseMatrix						*factMat;
	RowDenseMatrix						*oldFactMat;
	SolverTuckerLS<Tuple>				**solv;

	SecureDenseVector 					*coreTensor;
	double normTensor;
	

public:
	SolverTuckerFact(TensorIndexWithMultipleFile<Tuple>* inDbPtr)
		: ptrIndex(inDbPtr), factMat(nullptr), oldFactMat(nullptr), solv(nullptr), 
			coreTensor(nullptr), normTensor(-1)
	{
		init();
	}

	~SolverTuckerFact()
	{
		for (uint mode = 0; mode < Option::szMode; mode++)
			delete solv[mode];
		delete[] factMat;
		delete[] oldFactMat;
		delete[] solv;
	}

	void init()
	{
		factMat = new RowDenseMatrix[Option::szMode];
		oldFactMat = new RowDenseMatrix[Option::szMode];
		solv = new SolverTuckerLS<Tuple>*[Option::szMode];
		for (uint mode = 0; mode < Option::szMode; mode++){
			if (mode == 0)
				factMat[mode].init(Option::szDim[mode], Option::szRank[mode], false);
			else
				factMat[mode].init(Option::szDim[mode], Option::szRank[mode]);

			oldFactMat[mode].init(Option::szDim[mode], Option::szRank[mode]);

			if (M_SCAN == Option::flagMethod)
				solv[mode] = new SolverTransLS<Tuple>(mode, ptrIndex, factMat);

			else if (M_PLAIN == Option::flagMethod)
				solv[mode] = new SolverRowLS<Tuple>(mode, ptrIndex, factMat);
			else { 
				cout << "Hybrid version is not developed yet :) Sorry" << endl;
				exit(0);
				;
			}
			solv[mode]->init();
		}

	}

	void solve(){
		int it(0);

		double endTime;
		double newLoss(0), oldLoss, lossChange(0), eigValueChange(0);
		double modelLoss(0);
		double startTimeEval, evalTime(0);
		double startTime = GetTickCount();
		cout << "Learning starts with " << Option::szThread << " threads " << endl;
		if (TC_ITER == Option::flagCondTerminal && 0 == Option::maxIteration)
			return;
	
		cout << "Iteration\t Wall clock time\tloss change" << endl;
		bool loopFlag = true;
		while (loopFlag){
			//eigValueChange = 0;
			for (int mode = 0; mode < Option::szMode; mode++){
				solv[mode]->doIter();
				//eigValueChange += solv[mode]->doIter() / Option::szMode;
			}
			
			oldLoss = newLoss;
			startTimeEval = GetTickCount();
			if (L_CORE == Option::flagLoss){ // && TC_EPSILON == Option::flagCondTerminal){
				newLoss = getCoreLoss();
			} else if (L_MODEL == Option::flagLoss){
				newLoss = getModelLoss();
			}
			lossChange = abs(newLoss - oldLoss);
			evalTime += GetTickCount() - startTimeEval;

			it++;
			if (TC_ITER == Option::flagCondTerminal)
				loopFlag = it < Option::maxIteration;
			else if (TC_EPSILON == Option::flagCondTerminal)
				loopFlag = (lossChange > Option::epsilon) || it <= 1;
			else
				loopFlag = (it < Option::maxIteration && lossChange > Option::epsilon) || it <= 1;

			cout << "Iter#" << it << "\t\t" << (GetTickCount() - startTime)/1000 << "\t\t"
					<<  lossChange << "\t" << endl;
		}
		cout << "\nElapsed time summary - " << endl;
		cout << "\tdecomposition:\t" << (GetTickCount() - startTime)/1000 << endl;
		cout << "\tloss eval:\t" << evalTime/1000 << endl;
		return;
	}

	void run(){
		this->solve();
		this->storeModel(Option::modelFileName);
	}

	void storeModel(const char* modelFileName){

		fstream fp(modelFileName, fstream::out | fstream::trunc);

		fp << Option::szMode << " ";
		for (uint mode = 0; mode < Option::szMode; mode++)
			fp << Option::szDim[mode] << ":" << Option::szRank[mode] << " ";
		fp << endl;

		for (uint mode = 0; mode < Option::szMode; mode++){
			fp << mode << "th Factor matrix (" << Option::szDim[mode] << "-by-" << Option::szRank[mode] << endl;
			this->print(fp, mode);
		}

		fp.flush();
		fp.close();
	}

	void print(ostream& os, uint mode){
		os << factMat[mode];
	}

	double getModelLoss(){
		double lossSum = 0, diff = 0;
		const uint szWidth2 = 8, szWidth = 4;
		PRECISION *oldPtr, *newPtr;
		long uint i, szTotal;

		for (uint mit = 0; mit < Option::szMode; mit++){
			PRECISION sqauredLoss[szWidth];
			memset(sqauredLoss, 0, sizeof(PRECISION)*szWidth);

			oldPtr = oldFactMat[mit].ptrBulk;
			newPtr = factMat[mit].ptrBulk;
			szTotal = (long uint)Option::szDim[mit] * Option::szRank[mit];
			i = 0;

			#ifdef __AVX__
			__m256 xmmOldVec256, xmmNewVec256, xmmDiff256;
			__m256 xmmSum256 = _mm256_setzero_ps();
			for (i; i < szTotal - szWidth2 + 1; i = i + szWidth2){
				//_mm_prefetch(oldPtr + i + szWidth2, _MM_HINT_T0);
				//_mm_prefetch(newPtr + i + szWidth2, _MM_HINT_T0);
				//_mm_prefetch(oldPtr + i + szWidth2 + szWidth, _MM_HINT_T0);
				//_mm_prefetch(newPtr + i + szWidth2 + szWidth, _MM_HINT_T0);

				xmmOldVec256 = _mm256_load_ps(oldPtr + i);
				xmmNewVec256 = _mm256_load_ps(newPtr + i);
				xmmDiff256 = _mm256_sub_ps(xmmNewVec256, xmmOldVec256);
				xmmSum256 = _mm256_add_ps(xmmSum256, _mm256_mul_ps(xmmDiff256, xmmDiff256));
				_mm256_store_ps(oldPtr + i, xmmNewVec256);
			}

			__m128 *xmmPtr = (__m128*)&xmmSum256;
			__m128 xmmSum = _mm_add_ps(xmmPtr[0], xmmPtr[1]);

			if (i < szTotal - szWidth + 1){
				__m128 xmmOldVec, xmmNewVec, xmmDiff;
				xmmOldVec = _mm_load_ps(oldFactMat[mit].ptrBulk+i);
				xmmNewVec = _mm_load_ps(factMat[mit].ptrBulk+i);
				xmmDiff = _mm_sub_ps(xmmNewVec, xmmOldVec);
				xmmSum = _mm_add_ps(xmmSum, _mm_mul_ps(xmmDiff, xmmDiff));
				_mm_store_ps(oldFactMat[mit].ptrBulk+i, xmmNewVec);
				i = i + szWidth;
			}

			xmmSum = _mm_hadd_ps(xmmSum, xmmSum);
			xmmSum = _mm_hadd_ps(xmmSum, xmmSum);
			_mm_store_ps(sqauredLoss, xmmSum);
			#endif 

			for (i; i < szTotal; i++){
				diff = (oldPtr[i] - newPtr[i]);
				sqauredLoss[0] += diff*diff;
				oldFactMat[mit].ptrBulk[i] = factMat[mit].ptrBulk[i];
			}
	
			lossSum += sqrt(sqauredLoss[0]/Option::szRank[mit]);
		}

		return lossSum / Option::szMode;
	}

	double getCoreLoss(fstream* out = nullptr){

		uint lenCore = 1;
		SecureDenseVector *segCore = new SecureDenseVector[Option::szThread];
		double *segTensorLoss = nullptr;

		for (uint i = 0; i < Option::szMode; i++)
			lenCore = lenCore * Option::szRank[i];

		if (-1 == normTensor){
			segTensorLoss = new double[Option::szThread];
			memset(segTensorLoss, 0, sizeof(double)*Option::szThread);
			normTensor = 0;
		}

		Scanner::setMode(0);
		#pragma omp parallel for num_threads(Option::szThread)
		for (int tid = 0; tid < Option::szThread; tid++){
			segCore[tid].resize(lenCore, 0);
			Tuple temp;
			int lastRowID = -1;

			SecureDenseVector tempVec(lenCore);
			Scanner scan(ptrIndex);
			scan.rewind();
			while(scan.getNext(temp)){
				tempVec.resize(1, temp.value);
				if (nullptr != segTensorLoss)
					segTensorLoss[tid] += temp.value*temp.value;
				for (int mIter = Option::szMode; mIter > 0; mIter--)
					// Order preserving outer-product
					tempVec.outerProd(factMat[mIter-1].getVec(temp.dim[mIter-1]), Option::szRank[mIter-1], (out != nullptr));
				segCore[tid].acc(tempVec);
			}
		}

		// Aggregate losses
		if (nullptr!= segTensorLoss)
			normTensor = segTensorLoss[0];

		for (uint tid = 1; tid < Option::szThread; tid++){
			segCore[0].acc(segCore[tid]);
			if (nullptr!= segTensorLoss)
				normTensor += segTensorLoss[tid];
		}

		if (nullptr!= segTensorLoss){
			normTensor = sqrt(normTensor);
			delete[] segTensorLoss;
		}

		// Print core tensor 
		if (out != nullptr){
			(*out) << "CoreTensor" << endl;
			(*out) << segCore[0];
		}

		double coreSquare = innerProd(segCore[0], segCore[0]); 
		double coreLoss = 1 - sqrt(abs(normTensor*normTensor - coreSquare))/normTensor;

		//if (out == nullptr)
		//	cout << "( " << normTensor*normTensor << " - " << coreSquare << ") / " 
		//			<< normTensor << " = " << coreLoss << endl;

		delete[] segCore;

		if (std::isinf(coreLoss))
			return std::numeric_limits<PRECISION>::max();

		return coreLoss;
	}
};


#endif // #ifndef _ENGINE_TUCKER_FACT_
