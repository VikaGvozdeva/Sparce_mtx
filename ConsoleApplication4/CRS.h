#ifndef MODULES_COOMATRIX_INCLUDE_CRSMATRIX_H_
#define MODULES_COOMATRIX_INCLUDE_CRSMATRIX_H_
#include <inttypes.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>
using namespace std;

#define USE_DOUBLE
//#define USE_INT

#ifdef USE_INT
typedef int64_t INTTYPE;
#else
typedef int32_t INTTYPE;
#endif

#ifdef USE_DOUBLE
typedef double FPTYPE;
#define ZERO 0.0
#else
typedef float FPTYPE;
#define ZERO 0.0f
#endif

class CRSMatrix
{
public:
	INTTYPE N;
	INTTYPE NNZ;
	FPTYPE *val;
	INTTYPE *col_ind;
	INTTYPE *row_ptr;
	INTTYPE isSym = 0;

	CRSMatrix(INTTYPE _NNZ, INTTYPE _N, INTTYPE _isSym);
	CRSMatrix(const CRSMatrix& Matrix);
	~CRSMatrix();
	CRSMatrix& operator=(const CRSMatrix &Matrix);
	void ReadFromBinaryFile(char *filename);
	void WriteInBinaryFile(char* filename);
	void MatrixVectorMultCRS(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result);
	void SymMatrixVectorMultCRS(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result);
	friend ostream & operator<<(ostream &out, const CRSMatrix &Matrix);
	void DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N, INTTYPE _isSym);
	
};
#endif