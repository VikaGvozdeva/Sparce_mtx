#define _CRT_SECURE_NO_WARNINGS
#include "JD.h"
JDMatrix::JDMatrix(INTTYPE  _NNZ, INTTYPE _N, INTTYPE _MaxNNZ)
	{
//		int i;
		N = _N;
		NNZ = _NNZ;
		MaxNNZ = _MaxNNZ;
		if ((NNZ != 0) && (N != 0) && (MaxNNZ != 0))
		{
			jdiag = new FPTYPE[NNZ];
			perm = new INTTYPE[N];
			col_ind = new INTTYPE[NNZ];
			jd_ptr = new INTTYPE[MaxNNZ + 1];
			memset(jd_ptr, 0, (MaxNNZ + 1)  * sizeof(INTTYPE));
			memset(perm, 0, N * sizeof(INTTYPE));
			memset(jdiag, 0, NNZ * sizeof(FPTYPE));
			memset(col_ind, 0, NNZ * sizeof(INTTYPE));
		}
	}
	void JDMatrix::DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N, INTTYPE _MaxNNZ)
{
	if ((N != _N) || (NNZ != _NNZ) || (MaxNNZ != _MaxNNZ))
	{ 
		if ((NNZ != 0) && (N != 0) && (MaxNNZ != 0))
		{
			delete[] col_ind;
			delete[] jdiag;
			delete[] jd_ptr;
			delete[] perm;
		}

		if ((_N != 0) && (_NNZ != 0) && (_MaxNNZ != 0))
		{
			jdiag = new FPTYPE[_NNZ];
			perm = new INTTYPE[_N];
			col_ind = new INTTYPE[_NNZ];
			jd_ptr = new INTTYPE[_MaxNNZ + 1];
		}
		else
		{
			jdiag = 0;
			perm = 0;
			col_ind = 0;
			jd_ptr = 0;
		}
		memset(jd_ptr, 0, (MaxNNZ + 1)  * sizeof(INTTYPE));
		memset(perm, 0, N * sizeof(INTTYPE));
		memset(jdiag, 0, NNZ * sizeof(FPTYPE));
		memset(col_ind, 0, NNZ * sizeof(INTTYPE));
		N = _N;
		NNZ = _NNZ;
		MaxNNZ = _MaxNNZ;
	}
}
JDMatrix::~JDMatrix()
	{
		if ((N != 0) && (NNZ != 0))
			{
				delete[] perm;
				delete[] jdiag;
				delete[] col_ind;
				delete[] jd_ptr;
			}

	}
JDMatrix& JDMatrix::operator=(const JDMatrix &Matrix)
{
	if (this == &Matrix)
		return *this;
	DeleteAndAllocateMemory(Matrix.NNZ, Matrix.N, Matrix.MaxNNZ);
	memcpy(jd_ptr, Matrix.jd_ptr, sizeof(INTTYPE)*(MaxNNZ + 1));
	memcpy(perm, Matrix.perm, sizeof(INTTYPE)*N);
	memcpy(jdiag, Matrix.jdiag, sizeof(FPTYPE)*NNZ);
	memcpy(col_ind, Matrix.col_ind, sizeof(INTTYPE)*NNZ);
	return *this;
}
JDMatrix::JDMatrix(const JDMatrix &Matrix)
	{
		N = Matrix.N;
		NNZ = Matrix.NNZ;
		MaxNNZ = Matrix.MaxNNZ;
		jdiag = new FPTYPE[NNZ];
		perm = new INTTYPE[N];
		col_ind = new INTTYPE[NNZ];
		jd_ptr = new INTTYPE[MaxNNZ + 1];
		memcpy(jdiag, Matrix.jdiag, sizeof(FPTYPE)*NNZ);
		memcpy(perm, Matrix.perm, sizeof(INTTYPE)*N);
		memcpy(col_ind, Matrix.col_ind, sizeof(INTTYPE)*NNZ);
		memcpy(jd_ptr, Matrix.jd_ptr, sizeof(INTTYPE)*(MaxNNZ + 1));
	
	}
	void JDMatrix::ReadFromBinaryFile(char *filename)
	{
		FILE *JDmtx = NULL;
		int _N, _NNZ, _MaxNNZ;
		JDmtx = fopen(filename, "rb");
		if (JDmtx == NULL)
		{
			printf("Error opening file");
		}
		fread(&_N, sizeof(INTTYPE), 1, JDmtx);
		fread(&_NNZ, sizeof(INTTYPE), 1, JDmtx);
		fread(&_MaxNNZ, sizeof(INTTYPE), 1, JDmtx);
		DeleteAndAllocateMemory(_NNZ, _N, _MaxNNZ);
		fread(jdiag, sizeof(FPTYPE), NNZ, JDmtx);
		fread(col_ind, sizeof(INTTYPE), NNZ, JDmtx);
		fread(jd_ptr, sizeof(INTTYPE), MaxNNZ + 1, JDmtx);
		fread(perm, sizeof(INTTYPE), N, JDmtx);
		fclose(JDmtx);

	}

	void JDMatrix::WriteInBinaryFile(char* filename)
	{
		FILE *JDmtx = NULL;
		JDmtx = fopen(filename, "wb");
		if (JDmtx == NULL)
		{
			printf("Error opening file");
		}
		fwrite(&N, sizeof(INTTYPE), 1, JDmtx);
		fwrite(&NNZ, sizeof(INTTYPE), 1, JDmtx);
		fwrite(&MaxNNZ, sizeof(INTTYPE), 1, JDmtx);
		fwrite(jdiag, sizeof(FPTYPE), NNZ, JDmtx);
		fwrite(col_ind, sizeof(INTTYPE), NNZ, JDmtx);
		fwrite(jd_ptr, sizeof(INTTYPE), MaxNNZ + 1, JDmtx);
		fwrite(perm, sizeof(INTTYPE), N, JDmtx);
		fclose(JDmtx);
	}

	std::ostream & operator<<(ostream &out, const JDMatrix &Matrix)
	{
		out << "jdiag : " << endl;
		for (INTTYPE i = 0; i < Matrix.NNZ; i++)
		{
			out << Matrix.jdiag[i] << "  ";
		}
		out << endl;
		out << "col_ind : " << endl;
		for (INTTYPE i = 0; i < Matrix.NNZ; i++)
		{
			out << Matrix.col_ind[i] << "  ";
		}
		out << endl;
		out << "jd_ptr : " << endl;
		for (INTTYPE i = 0; i < Matrix.MaxNNZ + 1; i++)
		{
			out << Matrix.jd_ptr[i] << "  ";
		}
		out << endl;
		out << "perm : " << endl;
		for (INTTYPE i = 0; i < Matrix.N; i++)
		{
			out << Matrix.perm[i] << "  ";
		}
		out <<endl;
		return out;
	}

	void JDMatrix::MatrixVectorMultJD(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
	{
		int i = 0, j = 0, k = 0, NNZ = 0, tmp_ind = 0, maxval = 0, upper = 0;

		FPTYPE* temp = new FPTYPE[N];
		memset(result, 0, vec_N * sizeof(FPTYPE));

		int disp = 0;
//#pragma omp parallel for
		for (int j = 0; j < MaxNNZ; j++)
		{
//#pragma omp parallel for
			for (int i = 0; i < (jd_ptr[j + 1] - jd_ptr[j]); i++)
			{
				result[i] += jdiag[disp] * vec[col_ind[disp]];
				disp++;
			}
		}

		for (i = 0; i < N; i++)
		{
			temp[perm[i]] = result[i];
		}
		for (i = 0; i < N; i++)
		{
			result[i] = temp[i];
		}

		free(temp);
	}

	void JDMatrix::SymMatrixVectorMultJD(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
	{
		
	}

