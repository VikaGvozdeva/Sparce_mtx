#define _CRT_SECURE_NO_WARNINGS

#include "CRS.h"

CRSMatrix::CRSMatrix(INTTYPE _NNZ, INTTYPE _N, INTTYPE _isSym)
	{
		N = _N;
		NNZ = _NNZ;
		isSym = _isSym;
		if ((N != 0) && (NNZ != 0))
		{
			val = new FPTYPE[NNZ];
			col_ind = new INTTYPE[NNZ];
			row_ptr = new INTTYPE[N + 1];
			memset(val, 0, NNZ  * sizeof(FPTYPE));
			memset(col_ind, 0, N * sizeof(INTTYPE));
			memset(row_ptr, 0, (N+1) * sizeof(INTTYPE));
		}
		else
		{
			val = 0;
			col_ind = 0;
			row_ptr = 0;

		}
	}

void CRSMatrix::DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N, INTTYPE _isSym)
{
	if ((N != _N) || (NNZ != _NNZ))
	{ 
		if ((NNZ != 0) && (N != 0))
		{
			delete[] col_ind;
			delete[] val;
			delete[] row_ptr;
		}

		if ((_N != 0) && (_NNZ != 0))
		{
			val = new FPTYPE[_NNZ];
			col_ind = new INTTYPE[_NNZ];
			row_ptr = new INTTYPE[_N+1];
			memset(val, 0, _NNZ * sizeof(FPTYPE));
			memset(col_ind, 0, _N * sizeof(INTTYPE));
			memset(row_ptr, 0, (_N + 1) * sizeof(INTTYPE));
		}
		else
		{
			val = 0;
			col_ind = 0;
			row_ptr = 0;
		}
		
		N = _N;
		NNZ = _NNZ;
		isSym = _isSym;
	}
}
CRSMatrix::CRSMatrix(const CRSMatrix& Matrix)
	{
		N = Matrix.N;
		NNZ = Matrix.NNZ;
		val = new FPTYPE[NNZ];
		col_ind = new INTTYPE[NNZ];
		row_ptr = new INTTYPE[N + 1];
		memcpy(val, Matrix.val, sizeof(FPTYPE)*NNZ);
		memcpy(col_ind, Matrix.col_ind, sizeof(INTTYPE)*NNZ);
		memcpy(row_ptr, Matrix.row_ptr, sizeof(INTTYPE)*(N + 1));
	}

CRSMatrix::~CRSMatrix() 
{
	if ((N != 0) && (NNZ != 0))
	{
		delete[] col_ind;
		delete[] row_ptr;
		delete[] val;
	}
}

void CRSMatrix::ReadFromBinaryFile(char *filename)
	{
		FILE *CRSmtx = NULL;
		INTTYPE _N, _NNZ, _isSym;
		CRSmtx = fopen(filename, "rb");
		if (CRSmtx == NULL)
		{
			printf("Error opening file");
		}
		fread(&_N, sizeof(INTTYPE), 1, CRSmtx);
		fread(&_NNZ, sizeof(INTTYPE), 1, CRSmtx);
		fread(&_isSym, sizeof(INTTYPE), 1, CRSmtx);
		DeleteAndAllocateMemory(_NNZ, _N, _isSym);
		//CRSMatrix * Matrix = new CRSMatrix(NNZ, N);
		fread(val, sizeof(FPTYPE), NNZ, CRSmtx);
		fread(col_ind, sizeof(INTTYPE), NNZ, CRSmtx);
		fread(row_ptr, sizeof(INTTYPE), N + 1, CRSmtx);
		fclose(CRSmtx);
	}

	void CRSMatrix::WriteInBinaryFile(char* filename)
	{
		FILE *CRSmtx = NULL;
		CRSmtx = fopen(filename, "wb");
		if (CRSmtx == NULL)
		{
			printf("Error opening file");
		}
		fwrite(&N, sizeof(INTTYPE), 1, CRSmtx);
		fwrite(&NNZ, sizeof(INTTYPE), 1, CRSmtx);
		fwrite(&isSym, sizeof(INTTYPE), 1, CRSmtx);
		fwrite(val, sizeof(FPTYPE), NNZ, CRSmtx);
		fwrite(col_ind, sizeof(INTTYPE), NNZ, CRSmtx);
		fwrite(row_ptr, sizeof(INTTYPE), N + 1, CRSmtx);
		fclose(CRSmtx);
	}

	std::ostream & operator<<(ostream &out, const CRSMatrix &Matrix)
	{
		out << "Val : " << endl;
		for (INTTYPE i = 0; i < Matrix.NNZ; i++)
		{
			out << Matrix.val[i] << "  ";
		}
		out << endl;
		out << "col_ind : " << endl;
		for (INTTYPE i = 0; i < Matrix.NNZ; i++)
		{
			out << Matrix.col_ind[i] << "  ";
		}
		out << endl;
		out << "row_ptr : " << endl;
		for (INTTYPE i = 0; i < Matrix.N; i++)
		{
			out << Matrix.row_ptr[i] << "  ";
		}
		out << endl;
		return out;
	}

	CRSMatrix& CRSMatrix::operator=(const CRSMatrix &Matrix)
	{
	if (this == &Matrix)
		return *this;
	DeleteAndAllocateMemory(Matrix.NNZ, Matrix.N, Matrix.isSym);
	memcpy(val, Matrix.val, sizeof(FPTYPE)*NNZ);
	memcpy(col_ind, Matrix.col_ind, sizeof(INTTYPE)*NNZ);
	memcpy(row_ptr, Matrix.row_ptr, sizeof(INTTYPE)*(N + 1));
	return *this;

	}
	void CRSMatrix::MatrixVectorMultCRS(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
	{
		int i, j;
		double tmp;
		memset(result, 0, vec_N * sizeof(FPTYPE));

		if (vec_N == N)
		{
			#pragma omp parallel for
			for (i = 0; i < N; i++)
			{
				tmp = 0;
				//#pragma omp parallel for
				for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
				{
					tmp += val[j] * vec[col_ind[j]];
				}
				result[i] = tmp;
			}
		}
	}

	void CRSMatrix::SymMatrixVectorMultCRS(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
	{
		int i, j;
//		double tmp;
		int cur_row;
		memset(result, 0, vec_N * sizeof(FPTYPE));

		if (vec_N == N)
		{
			//#pragma omp parallel for
			for (i = 0; i < N; i++)
			{
				cur_row = i;
				//tmp = 0;
				//#pragma omp parallel for
				for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
				{
					if (cur_row == col_ind[j])
						result[i] += vec[cur_row] * val[j];
					else
					{
						result[i] += vec[col_ind[j]] * val[j];
						result[col_ind[j]] += vec[i] * val[j];
					}
				}
			}
		}
	}

