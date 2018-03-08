#define _CRT_SECURE_NO_WARNINGS

#include "CCS.h"
#include <omp.h>
CCSMatrix::CCSMatrix(INTTYPE _NNZ, INTTYPE _N)
	{
		N = _N;
		NNZ = _NNZ;
		if ((N != 0) && (NNZ != 0))
		{
			val = new FPTYPE[NNZ];
			row_ind = new INTTYPE[NNZ];
			col_ptr = new INTTYPE[N + 1];
			memset(val, 0, NNZ  * sizeof(FPTYPE));
			memset(row_ind, 0, N * sizeof(INTTYPE));
			memset(col_ptr, 0, (N+1) * sizeof(INTTYPE));
			}
		else
		{
			val = 0;
			row_ind = 0;
			col_ptr = 0;
		}
	}

void CCSMatrix::DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N)
{
	if ((N != _N) || (NNZ != _NNZ))
	{ 
		if ((NNZ != 0) && (N != 0))
		{
			delete[] row_ind;
			delete[] val;
			delete[] col_ptr;
		}

		if ((_N != 0) && (_NNZ != 0))
		{
			val = new FPTYPE[_NNZ];
			row_ind = new INTTYPE[_NNZ];
			col_ptr = new INTTYPE[_N+1];
		}
		else
		{
			val = 0;
			row_ind = 0;
			col_ptr = 0;
		}
		memset(val, 0, _NNZ  * sizeof(FPTYPE));
		memset(row_ind, 0, _N * sizeof(INTTYPE));
		memset(col_ptr, 0, (_N+1) * sizeof(INTTYPE));
		N = _N;
		NNZ = _NNZ;
	}
}

CCSMatrix::CCSMatrix(const CCSMatrix& Matrix)
	{
		N = Matrix.N;
		NNZ = Matrix.NNZ;
		val = new FPTYPE[NNZ];
		row_ind = new INTTYPE[NNZ];
		col_ptr = new INTTYPE[N + 1];
		memcpy(val, Matrix.val, sizeof(FPTYPE)*NNZ);
		memcpy(row_ind, Matrix.row_ind, sizeof(INTTYPE)*NNZ);
		memcpy(col_ptr, Matrix.col_ptr, sizeof(INTTYPE)*(N + 1));

	}

CCSMatrix::~CCSMatrix() 
{
		if ((N != 0) && (NNZ != 0))
		{
		delete[] row_ind;
		delete[] col_ptr;
		delete[] val;
	}
	}

	void CCSMatrix::ReadFromBinaryFile(char *filename)
	{
		FILE *CCSmtx = NULL;
		INTTYPE _N, _NNZ;
		CCSmtx = fopen(filename, "rb");
		if (CCSmtx == NULL)
		{
			printf("Error opening file");
		}
		fread(&_N, sizeof(INTTYPE), 1, CCSmtx);
		fread(&_NNZ, sizeof(INTTYPE), 1, CCSmtx);
		DeleteAndAllocateMemory(_NNZ, _N);
		fread(val, sizeof(FPTYPE), _NNZ, CCSmtx);
		fread(row_ind, sizeof(INTTYPE), _NNZ, CCSmtx);
		fread(col_ptr, sizeof(INTTYPE), _N + 1, CCSmtx);
		fclose(CCSmtx);
	}

	void CCSMatrix::WriteInBinaryFile(char* filename)
	{
		FILE *CCSmtx = NULL;
		CCSmtx = fopen(filename, "wb");
		if (CCSmtx == NULL)
		{
			printf("Error opening file");
		}
		fwrite(&N, sizeof(INTTYPE), 1, CCSmtx);
		fwrite(&NNZ, sizeof(INTTYPE), 1, CCSmtx);
		fwrite(val, sizeof(FPTYPE),NNZ, CCSmtx);
		fwrite(row_ind, sizeof(INTTYPE), NNZ, CCSmtx);
		fwrite(col_ptr, sizeof(INTTYPE), N + 1, CCSmtx);
		fclose(CCSmtx);
	}
	CCSMatrix& CCSMatrix::operator=(const CCSMatrix &Matrix)
	{
		if (this == &Matrix)
		return *this;
	DeleteAndAllocateMemory(Matrix.NNZ, Matrix.N);
	memcpy(val, Matrix.val, sizeof(FPTYPE)*NNZ);
	memcpy(row_ind, Matrix.row_ind, sizeof(INTTYPE)*NNZ);
	memcpy(col_ptr, Matrix.col_ptr, sizeof(INTTYPE)*(N + 1));
	return *this;

	}
	void CCSMatrix::MatrixVectorMultCCS(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
	{
		int i, j;
		bool flag=false;
		memset(result, 0, vec_N * sizeof(FPTYPE));
		if (vec_N == N)
		{	//omp_set_dynamic(0); 
			
			//static int a=1;
			#pragma omp parallel for
			for (i = 0; i < N; i++)
			{
		/*		if(flag==false)
			{
				//cout << "count of threads mult before -" << omp_get_num_threads() << endl;
				//omp_set_num_threads(16);
				cout << "count of threads mult -" << omp_get_num_threads() << endl;
				flag=true;
			}*/
			//cout << "count of threads mult in parallel-" << omp_get_num_threads() << endl;

				//#pragma omp parallel for
				for (j = col_ptr[i]; j < col_ptr[i + 1]; j++)
				{
					result[row_ind[j]] += vec[i] * val[j];
				}
			}
		}
	}

	void CCSMatrix::SymMatrixVectorMultCCS(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
	{
		INTTYPE i, j;
		//		double tmp;
		INTTYPE cur_row;
		memset(result, 0, vec_N * sizeof(FPTYPE));

		if (vec_N == N)
		{

			#pragma omp parallel for
			for (i = 0; i < N; i++)
			{
				cur_row = i;
				//tmp = 0;
				//#pragma omp parallel for
				for (j = col_ptr[i]; j < col_ptr[i + 1]; j++)
				{
					if (cur_row == row_ind[j])
						result[i] += vec[cur_row] * val[j];
					else
					{
						result[row_ind[j]] += vec[i] * val[j];
						result[i] += vec[row_ind[j]] * val[j];			
					}
				}
			}
		}
	}

	std::ostream & operator<<(ostream &out, const CCSMatrix &Matrix)
	{
		out << "Val : " << endl;
		for (INTTYPE i = 0; i < Matrix.NNZ; i++)
		{
			out << Matrix.val[i] << "  ";
		}
		out << endl;
		out << "row_ind : " << endl;
		for (INTTYPE i = 0; i < Matrix.NNZ; i++)
		{
			out << Matrix.row_ind[i] << "  ";
		}
		out << endl;
		out << "col_ptr : " << endl;
		for (INTTYPE i = 0; i < Matrix.N; i++)
		{
			out << Matrix.col_ptr[i] << "  ";
		}
		out << endl;
		return out;
	}
