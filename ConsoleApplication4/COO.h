#ifndef MODULES_COOMATRIX_INCLUDE_COOMATRIX_H_
#define MODULES_COOMATRIX_INCLUDE_COOMATRIX_H_
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <time.h>
#include <string.h>
#include <iostream>
#include <inttypes.h>
#include <cstring>
#include <vector>


using namespace std;

//#define USE_INT_64
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
#define ABS fabs
#else
typedef float FPTYPE;
#define ZERO 0.0f
#define ABS abs
#endif


#define MAX_LINE_LEN 1000
//int compare(const void * x1, const void * x2)
//{
//	return (*(int*)x1 - *(int*)x2);
//
//}

class CooMatrix  // ìàòðèöà â êîîðäèíàòíîì ôîðìàòå
{
public:

	INTTYPE N;
	INTTYPE NNZ;
	FPTYPE *val;
	INTTYPE *row_ind;
	INTTYPE *col_ind;
	INTTYPE *nnz_row;
	INTTYPE diag_cd;//number of different diagonals
	INTTYPE maxval_jd;//max NNZ
	INTTYPE diag_elem;//number of elements on the main diagonal
	INTTYPE isSym = 0;
	INTTYPE avg_nnz_per_row;
	INTTYPE max_nnz_per_row;
	INTTYPE min_nnz_per_row;
	FPTYPE sigma;

	CooMatrix(INTTYPE _NNZ, INTTYPE _N, INTTYPE _isSym); //êîíñòðóêòîð èíèöèàëèçàòîð
	CooMatrix(const CooMatrix &Matrix); //êîíñòðóêòîð êîïèðîâàíèÿ
	~CooMatrix();//äåñòðóêòîð
	CooMatrix& operator=(const CooMatrix &Matrix);   // ïðèñâàèâàíèå
	void ReadMatrix(char * filename); //÷òåíèå ìàòðèöû
	void ReadFromBinaryFile(char *filename); // ÷òåíèå èç áèíàðíîãî ôàéëà
	void WriteInBinaryFile(char *filename); //çàïèñü â áèíàðíûé ôàéë
	friend ostream & operator<<(ostream &out, const CooMatrix &Matrix);
	void DeleteAndAllocateMemory(INTTYPE _N, INTTYPE _NNZ, INTTYPE _isSym);
	void SortWithVectorsRow();
	void SortWithVectorsCol();
	void DiagCDMatrix(const CooMatrix& Matrix);
	void maxvalJDMatrix(const CooMatrix& Matrix);
	void calculate_params();
	//??
	//INTTYPE nonZeroRows(const COOMatrix& Matrix);
	void diagSLMatrix(const CooMatrix& Matrix);

};
#endif  //  COOMATRIX_H_