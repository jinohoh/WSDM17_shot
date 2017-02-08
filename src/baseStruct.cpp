#include "baseStruct.hpp"

#ifdef __SSE3__
#include <xmmintrin.h>
#include <mmintrin.h>
#include <pmmintrin.h>
#endif // #ifdef SSE

#ifdef __AVX__
#include <immintrin.h>
#endif // #ifdef AVX

/*
 * Records and DB
 */

uint &Record::szMode = Option::szMode;

Record::Record() : dim(nullptr){
	dim = new uint[szMode];
	memset(dim, 0, sizeof(uint)*szMode);
	value = 0;
}

Record::Record(const Record& temp){
	dim = new uint[szMode];
	memcpy(dim, temp.dim, sizeof(uint)*szMode);
	value = temp.value;
}

Record::~Record(){
	if (dim != nullptr)
		delete[] dim;
}

Record& Record::operator=(const Record& temp){
	memcpy(dim, temp.dim, sizeof(uint)*szMode);
	value = temp.value;
	return *this;
}

void Record::decode(char* iPtr){
	memcpy(dim, iPtr, sizeof(uint)*szMode);
	value = *((float*)(iPtr + sizeof(uint)*szMode));
}

void Record::encode(char* oPtr){
	memcpy(oPtr, dim, sizeof(uint)*szMode);
	(*((float*)(oPtr + sizeof(uint)*szMode))) = value;
}

void Record::print(){
	for (int i = 0; i < szMode; i++)
		cout << dim[i] << "\t";
	cout << value << endl;
}

/*
 * Vectors and in-memory structures
 */
PRECISION innerProd(const DenseVector& lhs, const DenseVector& rhs, const uint size){
	return innerProd(lhs.ptr, rhs.ptr, size);
}

PRECISION innerProd(const SecureDenseVector& lhs, const SecureDenseVector& rhs){
	assert(lhs.len == rhs.len);
	return innerProd(lhs.ptr, rhs.ptr, lhs.len); 
}

PRECISION	innerProd(const SecureDenseVector& lhs, PRECISION* rhs){
	return innerProd(lhs.ptr, rhs, lhs.len);
}

PRECISION innerProd(const PRECISION* lhs, const PRECISION* rhs, const uint size){

	// NOTE. the following codes are optimized for float 
	const int szWidth = 4;
	const int szWidth2 = szWidth * 2;
	long int i = 0;
	PRECISION sum[szWidth];
	memset(sum, 0, sizeof(PRECISION)*szWidth);

#ifdef __AVX__
	__m256 xmm256Lhs;
	__m256 xmm256Rhs;
	__m256 xmm256Sum = _mm256_setzero_ps();
	__m128* sumPtr;
	for (i = 0 ; i < (long int) size - szWidth2 + 1; i = i + szWidth2){
		//_mm_prefetch(lhs + i + szWidth2, _MM_HINT_T0);
		//_mm_prefetch(rhs + i + szWidth2, _MM_HINT_T0);
		//_mm_prefetch(lhs + i + szWidth2 + szWidth, _MM_HINT_T0);
		//_mm_prefetch(rhs + i + szWidth2 + szWidth, _MM_HINT_T0);
		xmm256Lhs = _mm256_load_ps(lhs + i);
		xmm256Rhs = _mm256_load_ps(rhs + i);
		xmm256Sum = _mm256_add_ps(xmm256Sum, _mm256_mul_ps(xmm256Lhs, xmm256Rhs));
	}
	sumPtr = (__m128*) &xmm256Sum;
	__m128 xmmSum = _mm_add_ps(sumPtr[0], sumPtr[1]);
	if ( i + szWidth < (int)size - 1) {
		__m128 xmmLhs;
		__m128 xmmRhs;
		xmmLhs = _mm_load_ps(lhs + i);
		xmmRhs = _mm_load_ps(rhs + i);
		xmmSum = _mm_add_ps(xmmSum, _mm_mul_ps(xmmLhs, xmmRhs));
		i = i + szWidth;
	}

	xmmSum = _mm_hadd_ps(xmmSum, xmmSum);
	xmmSum = _mm_hadd_ps(xmmSum, xmmSum);
	_mm_store_ps(sum, xmmSum);

	for (; i < size; i++)
		sum[0] += lhs[i] * rhs[i];

#else 
	for (i = 0; i < size; i++)
		sum[0] += lhs[i] * rhs[i];

#endif 

	return sum[0];
}

SecureDenseVector& SecureDenseVector::outerProd(DenseVector& lhs, const uint lenLhs, bool flagOrder){

	if (lenRsv < len * lenLhs){
		reserve(len * lenLhs);
	}

	PRECISION* tempPtr;
	long int j;

	if (!flagOrder){
		#ifdef __AVX__
		// NOTE. the following codes are optimized for float 
		const int szWidth2 = 8;
		const int szWidth = 4;
		uint cutLen = (len/szWidth2)*szWidth2;
		PRECISION remainder[szWidth2];
		memcpy(remainder, (ptr + cutLen), sizeof(PRECISION)*(len - cutLen));
	
		__m256 xmmScalar256;
		__m256 xmmVec256;
		__m128 xmmScalar128;
		__m128 xmmVec128;
	
		// AVX part
		if (cutLen > 0){
			for (uint i = 1; i < lenLhs; i++){
				tempPtr = ptr + (long uint)i*cutLen;
				if (0 == lhs[i]){
					memset(tempPtr, 0, sizeof(PRECISION)*cutLen);
				} else {
					xmmScalar256 = _mm256_set1_ps(lhs[i]);
					for (j = 0; j < cutLen; j = j + szWidth2){
						//_mm_prefetch(ptr + j + szWidth2, _MM_HINT_T0);
						//_mm_prefetch(ptr + j + szWidth2 + szWidth, _MM_HINT_T0);
						xmmVec256 = _mm256_load_ps(ptr + j);
						xmmVec256 = _mm256_mul_ps(xmmVec256, xmmScalar256);
						_mm256_store_ps(tempPtr + j, xmmVec256);
					}
				}
			}
			xmmScalar256 = _mm256_set1_ps(lhs[0]);
			for (j = 0; j < cutLen; j = j + szWidth2){
				xmmVec256 = _mm256_load_ps(ptr + j);
				_mm256_store_ps(ptr + j, _mm256_mul_ps(xmmVec256, xmmScalar256));
			}
		}
	
		// SSE part
		PRECISION* pivotPtr = ptr + (long uint)lenLhs * cutLen;
		cutLen = ((len - cutLen)/szWidth)*szWidth;
		if (cutLen > 0){
			for (uint i = 0; i < lenLhs; i++){
				tempPtr = pivotPtr + (long uint)i*cutLen;
				if (0 == lhs[i]){
					memset(tempPtr, 0, sizeof(PRECISION)*cutLen);
				} else {
					xmmScalar128 = _mm_set1_ps(lhs[i]);
					xmmVec128 = _mm_load_ps(remainder);
					_mm_store_ps(tempPtr, _mm_mul_ps(xmmVec128, xmmScalar128));
				}
			}
		}
	
		// plain part
		uint pivot = cutLen;
		pivotPtr = pivotPtr + (long uint)lenLhs * cutLen;
		uint szRemain = len % szWidth;
	
		if (szRemain != 0){
			for (uint i = 0; i < lenLhs; i++){
				tempPtr = pivotPtr + (long uint)i * szRemain;
				if (0 == lhs[i])
					memset(tempPtr, 0, sizeof(PRECISION)*szRemain);
				else {
					for (j = 0; j < szRemain; j++)
						tempPtr[j] = lhs[i] * remainder[pivot + j];
				}
			}
		}
		#else
		for (uint i = 1; i < lenLhs; i++){
			tempPtr = ptr + i*len;
			for (uint j = 0; j < len; j++)
			      tempPtr[j] = lhs[i] * ptr[j];
			
			for (uint j = 0; j < len; j++)
				ptr[i*len + j] = lhs[i] * ptr[j];
		}
		for (uint j = 0; j < len; j++)
			ptr[j] = lhs[0] * ptr[j];
		#endif
	}
	else {
		for (uint i = 1; i < lenLhs; i++){
			tempPtr = ptr + i*len;
			for (uint j = 0; j < len; j++)
			      tempPtr[j] = lhs[i] * ptr[j];
			
			for (uint j = 0; j < len; j++)
				ptr[i*len + j] = lhs[i] * ptr[j];
		}
		for (uint j = 0; j < len; j++)
			ptr[j] = lhs[0] * ptr[j];
	}
	
	len = len * lenLhs;
	return (*this);
}

void SecureDenseVector::reserve(uint length){
	if (lenRsv < length){
		lenRsv = length;
		PRECISION *oldptr = ptr;
		posix_memalign((void**)&ptr, ALIGN_WIDTH, (long long)sizeof(PRECISION)*lenRsv);
		// NOTE. len is always smaller than lenRsv and new length because resize does not decrease the space;
		memcpy(ptr, oldptr, len*sizeof(PRECISION));
		if (oldptr != nullptr)
			free(oldptr);
	}
}

void SecureDenseVector::resize(uint length){
	reserve(length);
	len = length;
}
void SecureDenseVector::resize(uint length, PRECISION value){
	resize(length);
	if (0 == value)
		memset(ptr, 0, sizeof(PRECISION)*len);
	else 
		for (uint i = 0; i < length; i++)
			ptr[i] = value;
}
void SecureDenseVector::initVal(PRECISION value){
	for (uint i = 0; i < len; i++)
		ptr[i] = value;
}

void SecureDenseVector::acc(PRECISION scalar, PRECISION* rhs){
	if (0 == scalar)
		return ;
	
	long int i = 0;
	const int szWidth2 = 8;
	const int szWidth = 4;

	#ifdef __AVX__
	__m256 xmmScalar256 = _mm256_set1_ps(scalar);
	__m256 xmmLhs256, xmmRhs256;
	__m128 xmmScalar = _mm_set1_ps(scalar);
	__m128 xmmLhs, xmmRhs;

	for (i = 0; i < (long int)len - szWidth2 + 1; i = i + szWidth2){
		//_mm_prefetch(ptr+ i + szWidth2, _MM_HINT_T0);
		//_mm_prefetch(rhs + i + szWidth2, _MM_HINT_T0);
		//_mm_prefetch(ptr + i + szWidth2 + szWidth, _MM_HINT_T0);
		//_mm_prefetch(rhs + i + szWidth2 + szWidth, _MM_HINT_T0);
		xmmLhs256 = _mm256_load_ps(ptr + i);
		xmmRhs256 = _mm256_load_ps(rhs + i);
		if (scalar != 1)
			xmmRhs256 = _mm256_mul_ps(xmmScalar256, xmmRhs256);
		xmmLhs256 = _mm256_add_ps(xmmLhs256, xmmRhs256);
		_mm256_store_ps(ptr + i, xmmLhs256);
	}

	if(i + szWidth < len - 1){
		xmmLhs = _mm_load_ps(ptr + i);
		xmmRhs = _mm_load_ps(rhs + i);
		if (scalar != 1)
			xmmRhs = _mm_mul_ps(xmmScalar, xmmRhs);
		xmmLhs = _mm_add_ps(xmmLhs, xmmRhs);
		_mm_store_ps(ptr + i, xmmLhs); 
 		i = i + szWidth;
	}
	#endif

	for (i; i < len; i++)
		ptr[i] += scalar * rhs[i];

	return;
}

void SecureDenseVector::acc(PRECISION scalar, SecureDenseVector& inVec){
	assert(len == inVec.len);
	acc(scalar, inVec.ptr);
}

void SecureDenseVector::acc(SecureDenseVector& inVec){
	PRECISION* rhs = inVec.ptr;
	assert(len == inVec.len);
	long int i = 0;

	#ifdef __AVX__
	const int szWidth2 = 8;
	const int szWidth = 4;

	__m256 xmmLhs256, xmmRhs256;
	__m128 xmmLhs, xmmRhs;

	for (i; i < (long int) len - szWidth2 + 1; i = i + szWidth2){
		xmmLhs256 = _mm256_load_ps(ptr + i);
		xmmRhs256 = _mm256_load_ps(rhs + i);
		xmmLhs256 = _mm256_add_ps(xmmLhs256, xmmRhs256);
		_mm256_stream_ps(ptr + i, xmmLhs256);
	}

	if (len - 1 > i + szWidth){
		xmmLhs = _mm_load_ps(ptr + i);
		xmmRhs = _mm_load_ps(rhs + i);
		xmmLhs = _mm_add_ps(xmmLhs, xmmRhs);
		_mm_stream_ps(ptr + i, xmmLhs);
		i = i + szWidth;
	}
	#endif

	for (i; i < len; i++)
		ptr[i] += rhs[i];

	return;
}

ostream& operator<< (ostream& os, SecureDenseVector& vec){

	if (vec.len > 0){
		os << "[ " << vec[0];
		for (uint i = 1; i < vec.len; i++){
			os << ", " << vec[i];
		}
		os <<" ] ";
	} else {
		os << "[] ";
	}

	return os;
}

void testInnerProd(){
	Vec a(10), b(10);
	for (int i = 0; i < 10; i++){
		a.ptr[i] = (9 - i);
		b.ptr[i] = i;
	}
	
	PRECISION innerRes = innerProd(a, b, 10);
	
	cout << a << endl;
	cout << b << endl;

	printf("Ideal value: %.4f,\tComputed value:%.4f\n", 120., innerRes);

	return;
}

void testOuterProd(){

	Vec a(10), b(12);
	for (int i = 0; i < 10; i++){
		a.ptr[i] = (9 - i);
	}
	for (int i = 0; i < 12; i++)
		b.ptr[i] = i;
	Vec d(2);
	d.ptr[0] = 1;
	d.ptr[1] = 0.1;
	
	cout << a << endl;
	cout << b << endl;
	Vec c(1);
	c.resize(1, 1);
	cout << c << endl;
	c.outerProd(a);
	cout << c << endl;
	c.outerProd(b);
	cout << c << endl;
	c.outerProd(d);
	cout << c << endl;
}

