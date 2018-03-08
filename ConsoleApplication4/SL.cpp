#define _CRT_SECURE_NO_WARNINGS
#include "SL.h"

SLMatrix::SLMatrix(INTTYPE  _NNZ, INTTYPE _N, INTTYPE _diag, INTTYPE _isSym)
	{
//		int i;
		N = _N;
		NNZ = _NNZ;
		diag = _diag;
		isSym = _isSym;
		if ((NNZ != 0) && (N != 0) &&(diag!=0))
		{
			adiag = new FPTYPE[N];
			iptr = new INTTYPE[N + 1];
			memset(adiag, 0, N * sizeof(FPTYPE));
			memset(iptr, 0, (N + 1) * sizeof(INTTYPE));
			if (isSym == 0)
			{
				autr = new FPTYPE[(NNZ - diag) / 2];
				altr = new FPTYPE[(NNZ - diag) / 2];
				jptr = new INTTYPE[(NNZ - diag) / 2];
				memset(autr, 0, ((NNZ - diag) / 2) * sizeof(FPTYPE));
				memset(altr, 0, ((NNZ - diag) / 2) * sizeof(FPTYPE));
				memset(jptr, 0, ((NNZ - diag) / 2) * sizeof(INTTYPE));
			}
			else if (isSym == 1)
			{
				altr = new FPTYPE[NNZ - diag];
				jptr = new INTTYPE[NNZ - diag];
				memset(altr, 0, (NNZ - diag) * sizeof(FPTYPE));
				memset(jptr, 0, (NNZ - diag) * sizeof(INTTYPE));
			}
		}
	}
SLMatrix::SLMatrix(const SLMatrix &Matrix)
{
		N = Matrix.N;
		NNZ = Matrix.NNZ;
		isSym = Matrix.isSym;
		adiag = new FPTYPE[N];
		iptr = new INTTYPE[N + 1];	
		memcpy(adiag, Matrix.adiag, sizeof(FPTYPE)*N);
		memcpy(iptr, Matrix.iptr, sizeof(INTTYPE)*(N + 1));
		
		if (isSym == 0)
		{
			autr = new FPTYPE[(NNZ - diag) / 2];
			altr = new FPTYPE[(NNZ - diag) / 2];
			jptr = new INTTYPE[(NNZ - diag) / 2];
			memcpy(autr, Matrix.autr, sizeof(FPTYPE)*((NNZ - diag) / 2));
			memcpy(altr, Matrix.altr, sizeof(FPTYPE)*((NNZ - diag) / 2));
			memcpy(jptr, Matrix.jptr, sizeof(INTTYPE)*((NNZ - diag) / 2));
		}
		else if (isSym == 1)
		{
			altr = new FPTYPE[NNZ - diag];
			jptr = new INTTYPE[NNZ - diag];
			memcpy(autr, Matrix.autr, sizeof(FPTYPE)*(NNZ - diag));
			memcpy(altr, Matrix.altr, sizeof(FPTYPE)*(NNZ - diag));
			memcpy(jptr, Matrix.jptr, sizeof(INTTYPE)*(NNZ - diag));
		}
}

SLMatrix& SLMatrix::operator=(const SLMatrix &Matrix)
{
		if (this == &Matrix)
		return *this;
	DeleteAndAllocateMemory(Matrix.NNZ, Matrix.N, Matrix.diag, Matrix.isSym);
	memcpy(adiag, Matrix.adiag, sizeof(FPTYPE)*N);
	memcpy(iptr, Matrix.iptr, sizeof(INTTYPE)*(N + 1));
	if (Matrix.isSym == 0)
	{
		memcpy(autr, Matrix.autr, sizeof(FPTYPE)*((NNZ - diag) / 2));
		memcpy(altr, Matrix.altr, sizeof(FPTYPE)*((NNZ - diag) / 2));
		memcpy(jptr, Matrix.jptr, sizeof(INTTYPE)*((NNZ - diag) / 2));
	}
	else if (Matrix.isSym == 1)
	{
		memcpy(altr, Matrix.altr, sizeof(FPTYPE)*(NNZ - diag));
		memcpy(jptr, Matrix.jptr, sizeof(INTTYPE)*(NNZ - diag));
	}
	
	return *this;
}
SLMatrix::~SLMatrix()
{
	if ((N != 0) && (NNZ != 0) && (diag != 0))
	{
		if (isSym == 0)
			delete[] autr;
		delete[] altr;
		delete[] iptr;
		delete[] jptr;
		delete[] adiag;
	}
}


void SLMatrix::ReadFromBinaryFile(char *filename)
	{

		FILE *SLmtx = NULL;
		int _N, _NNZ, _diag, _isSym;
		SLmtx = fopen(filename, "rb");
		if (SLmtx == NULL)
		{
			printf("Error opening file");
		}
		fread(&_N, sizeof(INTTYPE), 1, SLmtx);
		fread(&_NNZ, sizeof(INTTYPE), 1, SLmtx);
		fread(&_diag, sizeof(INTTYPE), 1, SLmtx);
		fread(&_isSym, sizeof(INTTYPE), 1, SLmtx);
		DeleteAndAllocateMemory(_NNZ, _N, _diag, _isSym);
		fread(adiag, sizeof(FPTYPE),N, SLmtx);
		fread(altr, sizeof(FPTYPE), (NNZ - diag) / 2, SLmtx);
		if (isSym ==1)
			fread(autr, sizeof(FPTYPE), (NNZ - diag) / 2, SLmtx);
		fread(jptr, sizeof(INTTYPE), (NNZ - diag) / 2, SLmtx);
		fread(iptr, sizeof(INTTYPE), N + 1, SLmtx);
		fclose(SLmtx);

	}
	void SLMatrix::WriteInBinaryFile(char* filename)
	{
		FILE *SLmtx = NULL;
		SLmtx = fopen(filename, "wb");
		if (SLmtx == NULL)
		{
			printf("Error opening file");
		}
		fwrite(&N, sizeof(INTTYPE), 1, SLmtx);
		fwrite(&NNZ, sizeof(INTTYPE), 1, SLmtx);
		fwrite(&diag, sizeof(INTTYPE), 1, SLmtx);
		fwrite(&isSym, sizeof(INTTYPE), 1, SLmtx);
		fwrite(adiag, sizeof(FPTYPE), N, SLmtx);
		fwrite(altr, sizeof(FPTYPE), (NNZ - diag) / 2, SLmtx);
		if (isSym == 0)
			fwrite(autr, sizeof(FPTYPE), (NNZ - diag) / 2, SLmtx);
		fwrite(jptr, sizeof(INTTYPE), (NNZ - diag) / 2, SLmtx);
		fwrite(iptr, sizeof(INTTYPE), N + 1, SLmtx);
		fclose(SLmtx);
	}

	void SLMatrix::MatrixVectorMultSL(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
	{
		int i, j;
		int m = 0;
		int up = 0, low = 0;
		memset(result, 0, vec_N * sizeof(FPTYPE));

		for (i = 0; i < N; i++)
		{
			result[i] = vec[i] * adiag[i];
		}
		if (isSym == 0)
		{
//			cout << "sym" << endl;
			//nonsymmetric algorithm
//#pragma omp parallel for
			for (i = 1; i < N; i++)
		{
//#pragma omp parallel for
			for (j = iptr[i]; j < iptr[i + 1]; j++)
			{
				result[i] += vec[jptr[up]] * altr[up];
				result[jptr[low]] += vec[i] * autr[low];
				up++; low++;
			}
		}
	}
		else if (isSym == 1)
		{
//#pragma omp parallel for
			for (i = 1; i < N; i++)
			{
//#pragma omp parallel for
					for (j = iptr[i]; j < iptr[i + 1]; j++)
					{
						result[i] += vec[jptr[up]] * altr[up];
						result[jptr[up]] += vec[i] * altr[up];
						up++; 
					}
			}
		}
	}

	std::ostream & operator<<(ostream &out, const SLMatrix &Matrix)
	{

		out << "adiag : " << endl;
		for (INTTYPE i = 0; i < Matrix.N; i++)
		{
			out << Matrix.adiag[i] << "  ";
		}
		out << endl;
		out << "altr : " << endl;
		for (INTTYPE i = 0; i < (Matrix.NNZ - Matrix.diag)/2 ; i++)
		{
			out << Matrix.altr[i] << "  ";
		}
		out << endl;
		if (Matrix.isSym == 0)
		{
			out << "autr : " << endl;
			for (INTTYPE i = 0; i < (Matrix.NNZ - Matrix.diag)/2; i++)
			{
				out << Matrix.autr[i] << "  ";
			}
			out << endl;
		}
		out << "jdptr : " << endl;
		for (INTTYPE i = 0; i < (Matrix.NNZ - Matrix.diag)/2; i++)
		{
			out << Matrix.jptr[i] << "  ";
		}
		out <<endl;
			out << "iptr : " << endl;
		for (INTTYPE i = 0; i < Matrix.N + 1; i++)
		{
			out << Matrix.iptr[i] << "  ";
		}
		out <<endl;
		return out;
	}
	void SLMatrix::DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N,  INTTYPE _diag, INTTYPE _isSym)
	{
		if ((N != _N) || (NNZ != _NNZ) || (diag != _diag))
	{ 
		if ((NNZ != 0) && (N != 0) && (diag != 0))
		{
			if (isSym == 0)
				delete[] autr;
			delete[] altr;
			delete[] iptr;
			delete[] jptr;
			delete[] adiag;
		}

		if ((_N != 0) && (_NNZ != 0) && (_diag != 0))
		{
			adiag = new FPTYPE[_N];
			iptr = new INTTYPE[_N + 1];
			memset(adiag, 0, N * sizeof(FPTYPE));
			memset(iptr, 0, (N + 1) * sizeof(INTTYPE));
			if (_isSym == 0)
			{
				autr = new FPTYPE[(_NNZ - _diag) / 2];
				altr = new FPTYPE[(_NNZ - _diag) / 2];
				jptr = new INTTYPE[(_NNZ - _diag) / 2];
				memset(autr, 0, ((NNZ - diag) / 2) * sizeof(FPTYPE));
				memset(altr, 0, ((NNZ - diag) / 2) * sizeof(FPTYPE));
				memset(jptr, 0, ((NNZ - diag) / 2) * sizeof(INTTYPE));
			}
			else if (_isSym == 1)
			{
				altr = new FPTYPE[_NNZ - _diag];
				jptr = new INTTYPE[_NNZ - _diag];
				memset(autr, 0, (NNZ - diag) * sizeof(FPTYPE));
				memset(altr, 0, (NNZ - diag) * sizeof(FPTYPE));
				memset(jptr, 0, (NNZ - diag) * sizeof(INTTYPE));
			}
		
		}
		else
		{
			adiag = 0;
			altr = 0;
			autr = 0;
			jptr = 0;
			iptr = 0;
		}
		//memset(adiag, 0, N * sizeof(FPTYPE));
		//memset(altr, 0, ((NNZ - diag) / 2) * sizeof(FPTYPE));
		//memset(autr, 0, ((NNZ - diag) / 2) * sizeof(FPTYPE));
		//memset(jptr, 0, ((NNZ - diag) / 2) * sizeof(INTTYPE));
		//memset(iptr, 0, (N + 1) * sizeof(INTTYPE));
		N = _N;
		NNZ = _NNZ;
		diag = _diag;
		isSym = _isSym;
	}
	}
