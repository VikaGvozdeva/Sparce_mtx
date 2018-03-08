#define _CRT_SECURE_NO_WARNINGS
//#pragma visibility 
#include "CD.h"
#include <iostream>
#include <new>
#include <climits>
#include <exception> 
int compare(const void * x1,const void * x2)
{
	return (*(int*)x1 - *(int*)x2);

}
CDMatrix::CDMatrix(INTTYPE  _NNZ, INTTYPE _N, INTTYPE _B)
{
	INTTYPE i;
	//	int j;
	N = _N;
	NNZ = _NNZ;
	B = _B;
	if ((N != 0) && (B != 0))
	{
		diag = new INTTYPE[B];
		val_data = new FPTYPE[B * N]; // n rows, B columns
	}
	else
	{
		val_data = 0;
	}
	if (N != 0)
	{
		val = new FPTYPE*[N]; // n pointers to rows
		for (int i = 0; i < N; i++)
			val[i] = val_data + i * B;
	}
	else
	{
		val = 0;
	}

}

CDMatrix::~CDMatrix()
{
	if ((B != 0) && (N != 0))
	{
		delete[] val_data;
		delete[] diag;
	}
	if (N != 0)
	{
		delete[] val;
	}
}
CDMatrix::CDMatrix(const CDMatrix &Matrix)
{
	N = Matrix.N;
	NNZ = Matrix.NNZ;
	B = Matrix.B;
	diag = new INTTYPE[B];
	val_data = new FPTYPE[B * N];
	memcpy(val_data, Matrix.val_data, sizeof(FPTYPE)* B * N);
	memcpy(diag, Matrix.diag, sizeof(INTTYPE)*B);
	val = new FPTYPE*[N];
	memcpy(val, Matrix.val, sizeof(FPTYPE*)*N);
}

CDMatrix& CDMatrix::operator=(const CDMatrix &Matrix)
{
	if (this == &Matrix)
		return *this;
	DeleteAndAllocateMemory(Matrix.NNZ, Matrix.N, Matrix.B);
	memcpy(diag, Matrix.diag, sizeof(INTTYPE)*B);
	memcpy(val_data, Matrix.val_data, sizeof(FPTYPE)*N*B);
	memcpy(val, Matrix.val, sizeof(FPTYPE*)*N);
	return *this;
}
void CDMatrix::DeleteAndAllocateMemory(INTTYPE _NNZ, INTTYPE _N, INTTYPE _B)
{

	if ((N != _N) || (B != _B))
	{
		if ((B != 0) && (N != 0))
		{
			delete[] val_data;
			delete[] diag;
			delete[] val;
		}

		if ((_N != 0) && (_B != 0))
		{

cout<<"Memory:"<<endl;
cout<<"N :"<<_N<<endl;
cout<<"B :"<<_B<<endl;

 try
  {
    val_data = new FPTYPE[_N*_B];
  }
  catch (std::bad_alloc& ba)
  {
    std::cerr << "bad_alloc caught: " << ba.what() << '\n';
  }

/*			try
			{
				val_data = new FPTYPE[_N*_B];
			}
			catch(const std::bad_array_new_length &e) {
        cout <<"Memory problem - "<< e.what() << endl;
			//val_data = new FPTYPE[_N*_B];
		}*/
	}
		else
		{
			val_data = 0;
		}
		if (_N != 0)
		{
			val = new FPTYPE*[_N];
			for (INTTYPE i = 0; i < _N; i++)
			{
				val[i] = val_data + i*_B;
			}
			diag = new INTTYPE[_B];
		}
		else
		{
			val = 0;
		}
		//for (INTTYPE i = 0; i < _N; i++)
		//{
		//	val[i] = val_data + i*_B;
		//}
	}
	N = _N;
	NNZ = _NNZ;
	B = _B;
}
void CDMatrix::ReadFromBinaryFile(char *filename)
{
	FILE *CDmtx = NULL;
	INTTYPE _N, _NNZ, _B;
	CDmtx = fopen(filename, "rb");
	if (CDmtx == NULL)
	{
		printf("Error opening file");
	}
	fread(&_N, sizeof(INTTYPE), 1, CDmtx);
	fread(&_NNZ, sizeof(INTTYPE), 1, CDmtx);
	fread(&_B, sizeof(INTTYPE), 1, CDmtx);
	DeleteAndAllocateMemory(_NNZ, _N, _B);
	fread(diag, sizeof(INTTYPE), B, CDmtx);
	fread(val_data, sizeof(FPTYPE), N*B, CDmtx);
	for (INTTYPE i = 0; i < N; i++)
	{
		val[i] = val_data + i*B;
	}
	fclose(CDmtx);

}
void CDMatrix::FillDiagArray(const CooMatrix &Matrix)
{
	bool flag = false;
	INTTYPE p = 0, tmp_ind;
	//vector<INTTYPE> diag_ind(Matrix.diag_cd, (0 - Matrix.N - 1));
	memset(diag, 0, B * sizeof(INTTYPE));
	vector<INTTYPE> diag_ind;
	vector<INTTYPE>::iterator it;
	tmp_ind = Matrix.col_ind[0] - Matrix.row_ind[0];
	diag_ind.push_back(tmp_ind);
	//p++;
	//#pragma omp parallel for
	for (int i = 1; i < Matrix.NNZ; i++)
	{
		tmp_ind = Matrix.col_ind[i] - Matrix.row_ind[i];
		// bool falg = false;
		// for (int j=0; j<diag_ind.size();j++)
		// {
		// 	if (diag_ind[j]==tmp_ind)
		// 		flag=true;
		// }
		// if (flag==false)
		// 	diag_ind.push_back(tmp_ind);

		if (find(diag_ind.begin(), diag_ind.end(), tmp_ind) == diag_ind.end())
			diag_ind.push_back(tmp_ind);
	}
	for (int i = 0; i < diag_ind.size(); i++)
	{
		diag[i] = diag_ind[i];
	}

	qsort(diag, diag_ind.size(), sizeof(INTTYPE), compare);
	//	for (int i =0; i < diag_ind.size();i++)
	cout <<"diag vector - " <<diag_ind.size() <<endl;
	cout <<"diag B - " <<B <<endl;
}


void CDMatrix::WriteInBinaryFile(char* filename)
{
	FILE *CDmtx = NULL;
	CDmtx = fopen(filename, "wb");
	if (CDmtx == NULL)
	{
		printf("Error opening file");
	}
	fwrite(&N, sizeof(INTTYPE), 1, CDmtx);
	fwrite(&NNZ, sizeof(INTTYPE), 1, CDmtx);
	fwrite(&B, sizeof(INTTYPE), 1, CDmtx);
	fwrite(diag, sizeof(INTTYPE), B, CDmtx);
	fwrite(val_data, sizeof(FPTYPE), N * B, CDmtx);
	fclose(CDmtx);

}
void CDMatrix::MatrixVectorMultCD(FPTYPE *vec, INTTYPE vec_N, FPTYPE *result)
{
	INTTYPE i, j;
	INTTYPE tmp;
	INTTYPE max;

	//int B = Matrix->B;
	memset(result, 0, sizeof(FPTYPE)*N);
#pragma omp parallel for
	for (j = 0; j < B; j++) {
		tmp = diag[j];
		if (0 >(0 - tmp))
			i = 0;
		else i = 0 - tmp;

		if (0 > tmp)
			max = N - 0;
		else max = N - tmp;
		//#pragma omp parallel for
		for (i; i < max; i++)
			result[i] += val[i][j] * vec[tmp + i];

	}
}

void CDMatrix::SymMatrixVectorMultCD(FPTYPE * vec, INTTYPE vec_N, FPTYPE * result)
{
	INTTYPE up, low, j;
	INTTYPE tmp_low, tmp_up;
	INTTYPE max;
	INTTYPE k;
	bool main_diag = true;
	INTTYPE* diag_up = new INTTYPE[B];
	memcpy(diag_up, diag, sizeof(INTTYPE)*B);
	for (int i = 0; i < B; i++)
	{
		diag_up[i] *= -1;
	}
	if (diag_up[B - 1] != 0)
	{
		main_diag = false;
	}
	qsort(diag_up, B, sizeof(INTTYPE), compare);
	memset(result, 0, sizeof(FPTYPE)*N);

	//#pragma omp parallel for
	for (j = 0; j < B; j++) {
		tmp_low = diag[j];
		low = 0 - tmp_low;
		max = N;
		//#pragma omp parallel for
		for (low; low < max; low++)
			result[low] += val[low][j] * vec[tmp_low + low];
	}
	if (main_diag == false)
		j = 0;
	else
		j = 1;

	//#pragma omp parallel for
	for (j; j < B; j++)
	{
		tmp_up = diag_up[j];
		up = 0;
		max = N - tmp_up;
		//#pragma omp parallel for
		for (up; up < max; up++)
			result[up] += val[up][j] * vec[tmp_up + up];
	}
	delete[]  diag_up;
}

std::ostream & operator<<(ostream &out, const CDMatrix &Matrix)
{
	//print CD matrix
	return out;
}