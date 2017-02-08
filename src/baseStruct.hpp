#ifndef _BASESTRUCT_
#define _BASESTRUCT_
#include "Macros.hpp"
#include "Option.hpp"
#include <cstring> // memcpy

using namespace std;

/*
 * Primitives for File accessing
 */

struct RecordID {
	uint pageNO;
	uint slotNO;
};

struct Record {
	uint		*dim;
	PRECISION	value;
	static uint &szMode;
	static uint getSize(){ return szMode*sizeof(uint) + sizeof(PRECISION); }
	Record();
	Record(const Record&);
	~Record();
	Record& operator=(const Record&);
	void encode(char* oPtr);
	void decode(char* iPtr);
	void print();
};


/***************************************************************
 * Plain Vector
 **************************************************************/
struct DenseVector {
	PRECISION *ptr;
	DenseVector(): ptr(nullptr) {}
	~DenseVector() {
		if (ptr != nullptr)
			free(ptr);
	}
	PRECISION& operator[](uint index){ return ptr[index]; }
};
typedef DenseVector DenseVector;


/***************************************************************
 * Secure Vector
 **************************************************************/
struct SecureDenseVector : public DenseVector {
	uint lenRsv;
	uint len;

	SecureDenseVector(): DenseVector(), lenRsv(0), len(0) {}
	SecureDenseVector(uint length): len(length), lenRsv(length) {
		if (length > 0)
			posix_memalign((void**)&ptr, ALIGN_WIDTH, (long long)sizeof(PRECISION)*lenRsv);
		memset(ptr, 0, sizeof(PRECISION)*len);
	}

	void reserve(uint length);
	void resize(uint length);
	void resize(uint length, PRECISION value);
	void initVal(PRECISION value);

	void acc(PRECISION scalar, PRECISION* inPtr);
	void acc(PRECISION scalar, SecureDenseVector& inVec);
	void acc(SecureDenseVector& inVec);

	SecureDenseVector& outerProd(DenseVector& lhs, const uint lenLhs, bool flagOrder = false);
	SecureDenseVector& outerProd(SecureDenseVector& lhs){
		return outerProd((DenseVector&) lhs, lhs.len); 
	}
};

typedef SecureDenseVector Vec;

ostream& operator<<(ostream& os, SecureDenseVector&);

PRECISION 	innerProd(const DenseVector&, const DenseVector&, const uint size);
PRECISION	innerProd(const SecureDenseVector&, const SecureDenseVector&);
PRECISION	innerProd(const SecureDenseVector&, PRECISION*);
PRECISION	innerProd(const PRECISION*, const PRECISION*, const uint size);


/*
 * Test codes
 */

void testInnerProd();
void testOuterProd();

#endif // #ifndef _BASESTRUCT_
