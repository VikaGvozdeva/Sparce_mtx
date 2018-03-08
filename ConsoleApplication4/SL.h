#ifndef MODULES_COOMATRIX_INCLUDE_SLMATRIX_H_
#define MODULES_COOMATRIX_INCLUDE_SLMATRIX_H_
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include <string.h>
#include <iostream>
#include <inttypes.h>
#include "COO.h"

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

//#define MAX_LINE_LEN 1000000
class SLMatrix
{
public:

	INTTYPE N;
	INTTYPE NNZ;
	INTTYPE diag;
	FPTYPE *adiag;
	FPTYPE *altr;
	FPTYPE *autr;
	INTTYPE *jptr;
	INTTYPE *iptr;
	INTTYPE isSym;


	SLMatrix(INTTYPE  _NNZ, INTTYPE _N, INTTYPE _diag, INTTYPE _isSym);
	SLMatrix(const SLMatrix &Matrix); 
	SLMatrix& operator=(const SLMatrix &Matrix);
	~SLMatrix();
	void ReadFromBinaryFile(char *filename);
	void WriteInBinaryFile( char* filename);
	void MatrixVectorMultSL(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result);
	friend ostream & operator<<(ostream &out, const SLMatrix &Matrix);
	void DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N,  INTTYPE _diag, INTTYPE _isSym);

};
#endif