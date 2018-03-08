#pragma once
#ifndef MODULES_COOMATRIX_INCLUDE_JDMATRIX_H_
#define MODULES_COOMATRIX_INCLUDE_JDMATRIX_H_
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include <string.h>
#include <iostream>
#include <inttypes.h>

using namespace std;
//#define USE_INT
#define USE_DOUBLE

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
class JDMatrix
{
public:
	INTTYPE NNZ;
	INTTYPE N;
	INTTYPE MaxNNZ;
	INTTYPE *perm;
	FPTYPE *jdiag;
	INTTYPE *col_ind;
	INTTYPE *jd_ptr;
	INTTYPE isSym = 0;

	JDMatrix(INTTYPE  _NNZ, INTTYPE _N, INTTYPE _MaxNNZ);	
	JDMatrix(const JDMatrix &Matrix);
	~JDMatrix();
	void ReadFromBinaryFile(char *filename);	
	void WriteInBinaryFile(char *filename);		
	void MatrixVectorMultJD(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result);
	void SymMatrixVectorMultJD(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result);
	JDMatrix& operator=(const JDMatrix &Matrix);
	friend ostream & operator<<(ostream &out, const JDMatrix &Matrix);
	void DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N,  INTTYPE _MaxNNZ);
	
};
#endif
