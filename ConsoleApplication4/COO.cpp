#define _CRT_SECURE_NO_WARNINGS
#include "COO.h"
#include "CD.h"

CooMatrix::CooMatrix(INTTYPE _NNZ, INTTYPE _N, INTTYPE _isSym)
{
	N = _N;
	NNZ = _NNZ;
	isSym = _isSym;


	if (NNZ != 0)
	{

cout<<"NNZ coo initialize "<<NNZ<<endl;
try {
		val = new FPTYPE[NNZ];
		col_ind = new INTTYPE[NNZ];
		row_ind = new INTTYPE[NNZ];
	}
  catch (std::bad_alloc& ba)
  {
    std::cerr << "coo init bad_alloc caught: " << ba.what() << '\n';
  }

	}
	if (N != 0)
	{
		nnz_row = new INTTYPE[N];
		memset(nnz_row, 0, sizeof(INTTYPE)*N);
	}
}

CooMatrix::CooMatrix(const CooMatrix &Matrix) 
{
	N = Matrix.N;
	NNZ = Matrix.NNZ;
	isSym = Matrix.isSym;
	val = new FPTYPE[NNZ];
	col_ind = new INTTYPE[NNZ];
	row_ind = new INTTYPE[NNZ];
	nnz_row = new INTTYPE[N];
	memcpy(nnz_row, Matrix.nnz_row, sizeof(INTTYPE)*N);
	memcpy(val, Matrix.val, sizeof(FPTYPE)*NNZ);
	memcpy(col_ind, Matrix.col_ind, sizeof(INTTYPE)*NNZ);
	memcpy(row_ind, Matrix.row_ind, sizeof(INTTYPE)*NNZ);
}

CooMatrix :: ~CooMatrix()
{
	if (N != 0)
		delete[] nnz_row;
	if (NNZ != 0)
	{
		delete[] val;
		delete[] col_ind;
		delete[] row_ind;
	}
}

CooMatrix& CooMatrix :: operator=(const CooMatrix &Matrix)    
{
	if (this == &Matrix)
		return *this;
	DeleteAndAllocateMemory(Matrix.N, Matrix.NNZ, Matrix.isSym);
	memcpy(nnz_row, Matrix.nnz_row, sizeof(INTTYPE)*N);
	memcpy(val, Matrix.val, sizeof(FPTYPE)*NNZ);
	memcpy(col_ind, Matrix.col_ind, sizeof(INTTYPE)*NNZ);
	memcpy(row_ind, Matrix.row_ind, sizeof(INTTYPE)*NNZ);
	return *this;
}

void CooMatrix::ReadMatrix(char * filename) 
{
	FILE* f;
	char* line, *line_;
	char* p = NULL;
	char* sym = NULL;
	INTTYPE NNZ_file=0, diag_elem=0;
	INTTYPE _N, _NNZ, _isSym;
	f = fopen(filename, "r");
	if (f == NULL)
		printf("%s- File Not Found!\n", filename);
	//line_ = new char[200];
	//int i = 0;
	//do
	//	line_[i++] = fgetc(f);
	//while
	//	(line_[i] != '-');

	////sym = strtok(line_, "symmetric");
	//sym = strstr(line_, "symmetric");
	//if (sym != NULL) _isSym = 1;

	line = new char[MAX_LINE_LEN];

	char* line_sym = new char[MAX_LINE_LEN];
	fgets(line_sym, MAX_LINE_LEN, f);
	if (line_sym[0] == '%')
	{
		if (strstr(line_sym, "symmetric")!=NULL)
		{
			_isSym = 1;
		}
		else _isSym = 0;
	}

	do
		fgets(line, MAX_LINE_LEN, f);
	while (line[0] == '%');

	p = strtok(line, " ");
	_N = atoi(p);
	p = strtok(NULL, " ");
	p = strtok(NULL, " ");
	_NNZ = atoi(p);
	cout<<"before read elements"<<endl;
	//DeleteAndAllocateMemory(_N, _NNZ, _isSym);
	INTTYPE* _nnz_row = new INTTYPE[_N];
	memset(_nnz_row,0,sizeof(INTTYPE)*_N);
	vector <INTTYPE> _row_ind;
	vector <INTTYPE> _col_ind;
	vector <FPTYPE> _val;
	cout<<"vectors created"<<endl;

	INTTYPE temp_col_ind, temp_row_ind;
	FPTYPE temp_val;
	if (_isSym == 1)
	{
		for (INTTYPE i = 0; i < _NNZ; i++)
		{
			fgets(line, MAX_LINE_LEN, f);
			p = strtok(line, " ");
			//_row_ind.push_back(atoi(p) - 1);
			temp_row_ind = atoi(p) - 1;
			p = strtok(NULL, " ");
			//_col_ind.push_back(atoi(p) - 1);
			temp_col_ind = atoi(p) - 1;
			p = strtok(NULL, " ");
			temp_val = atof(p);
			if ((temp_col_ind == temp_row_ind) && (temp_val != ZERO))
			{
				_val.push_back(temp_val);
				_row_ind.push_back(temp_row_ind);
				_col_ind.push_back(temp_col_ind);
				_nnz_row[temp_row_ind]++;
				NNZ_file++;
				diag_elem++;
			}
			else if ((temp_col_ind != temp_row_ind) && (temp_val != ZERO))
			{
				_val.push_back(temp_val);
				_val.push_back(temp_val);
				_row_ind.push_back(temp_row_ind);
				_row_ind.push_back(temp_col_ind);
				_col_ind.push_back(temp_col_ind);
				_col_ind.push_back(temp_row_ind);
				_nnz_row[temp_col_ind]++;
				_nnz_row[temp_row_ind]++;
				NNZ_file += 2;
			}
		}
		delete[]line;
		delete[]line_sym;
		fclose(f);
		cout<<"elements were red"<<endl;

		DeleteAndAllocateMemory(_N, NNZ_file, _isSym);
		memcpy(nnz_row, _nnz_row, sizeof(INTTYPE)*N);
		for (INTTYPE j = 0; j < NNZ_file; j++)
		{
			val[j] = _val[j];
			col_ind[j] = _col_ind[j];
			row_ind[j] = _row_ind[j];
		}
	}
	else if (isSym==0)
	{
		DeleteAndAllocateMemory(_N, _NNZ, _isSym);
		for (INTTYPE i = 0; i < _NNZ; i++)
		{
			fgets(line, MAX_LINE_LEN, f);
			p = strtok(line, " ");
			row_ind[i] = atoi(p) - 1;
			nnz_row[row_ind[i]]++;
			p = strtok(NULL, " ");
			col_ind[i] = atoi(p) - 1;
			p = strtok(NULL, " ");
			val[i] = atof(p);
		}
		delete[]line;
		fclose(f);
	}
	delete[] _nnz_row;
}

void CooMatrix::ReadFromBinaryFile(char *filename)  
{
	FILE *COOmtx = NULL;
	INTTYPE _N, _NNZ, _isSym;
	COOmtx = fopen(filename, "rb");
	if (COOmtx == NULL)
	{
		printf("Error opening file");
	}
	fread(&_N, sizeof(INTTYPE), 1, COOmtx);
	fread(&_NNZ, sizeof(INTTYPE), 1, COOmtx);
	fread(&_isSym, sizeof(INTTYPE), 1, COOmtx);
	DeleteAndAllocateMemory(_N, _NNZ, _isSym);
	fread(val, sizeof(FPTYPE), NNZ, COOmtx);
	fread(col_ind, sizeof(INTTYPE), NNZ, COOmtx);
	fread(row_ind, sizeof(INTTYPE), NNZ, COOmtx);
	fread(nnz_row, sizeof(INTTYPE), N, COOmtx);
	fclose(COOmtx);
}

void CooMatrix::WriteInBinaryFile(char *filename)  
{
	FILE *COOmtx = NULL;
	COOmtx = fopen(filename, "wb");
	if (COOmtx == NULL)
	{
		printf("Error opening file");
	}
	fwrite(&N, sizeof(INTTYPE), 1, COOmtx);
	fwrite(&NNZ, sizeof(INTTYPE), 1, COOmtx);
	fwrite(&isSym, sizeof(INTTYPE), 1, COOmtx);
	fwrite(val, sizeof(FPTYPE), NNZ, COOmtx);
	fwrite(col_ind, sizeof(INTTYPE), NNZ, COOmtx);
	fwrite(row_ind, sizeof(INTTYPE), NNZ, COOmtx);
	fwrite(nnz_row, sizeof(INTTYPE), N, COOmtx);
	fclose(COOmtx);
}

std::ostream & operator<<(ostream &out, const CooMatrix &Matrix)
{
	for (INTTYPE i = 0; i < Matrix.NNZ; i++)
	{
		out << "(" << Matrix.val[i] << "," << Matrix.row_ind[i] << "," << Matrix.col_ind[i] << ") ,";
	}
	out << endl;
	return out;
}
void CooMatrix::DeleteAndAllocateMemory(INTTYPE _N, INTTYPE _NNZ, INTTYPE _isSym)
{
	if (N != _N)
	{
		if (N != 0)
			delete[]nnz_row;
		nnz_row = new INTTYPE[_N];
		memset(nnz_row, 0, sizeof(INTTYPE)*_N);
		N = _N;
	}
	if (NNZ != _NNZ)
	{
		if (NNZ != 0)
		{
			cout<<"NNZ coo initialize "<<NNZ<<endl;
try {
		val = new FPTYPE[NNZ];
		col_ind = new INTTYPE[NNZ];
		row_ind = new INTTYPE[NNZ];
	}
  catch (std::bad_alloc& ba)
  {
    std::cerr << "coo allocate bad_alloc caught: " << ba.what() << '\n';
  }
			delete[]val;
			delete[]col_ind;
			delete[]row_ind;
		}
		val = new FPTYPE[_NNZ];
		col_ind = new INTTYPE[_NNZ];
		row_ind = new INTTYPE[_NNZ];
		NNZ = _NNZ;
		isSym = _isSym;
	}
}


	void CooMatrix::maxvalJDMatrix(const CooMatrix& Matrix)
	{
		maxval_jd = Matrix.nnz_row[0];
		for (int i = 1; i < N; i++)
		{
			if (Matrix.nnz_row[i] > maxval_jd)
			{
				maxval_jd = Matrix.nnz_row[i];
			}
		}
	}

	void CooMatrix::diagSLMatrix(const CooMatrix& Matrix)
	{
		diag_elem =  0;
		for (int i = 0; i < NNZ; i++)
		{
			if (Matrix.row_ind[i] == Matrix.col_ind[i])
				diag_elem++;
		}
	}

	void CooMatrix::DiagCDMatrix(const CooMatrix& Matrix)
	{
		int i = 0, j = 0, p = 0, k = 0, l = 0, NNZ = 0, N = 0, diag_ind = 0, B = 0, tmp_ind = 0, m = 0, diag_numb = 0;
		NNZ = Matrix.NNZ;
		N = Matrix.N;

		bool flag;

		INTTYPE *temp = new INTTYPE[2 * N - 1];

		for (int i = 0; i < 2 * N - 1; i++)
		{
			//diag[i] = N + 1;
			temp[i] = N + 1;
		}
		for (int i = 0; i < NNZ; i++)
		{
			tmp_ind = Matrix.col_ind[i] - Matrix.row_ind[i];

			for (int j = 0; j < 2 * N - 1; j++)
			{
				//if ((tmp_ind) != diag[j])
				if ((tmp_ind) != temp[j])
				{
					flag = true;
				}
				else
				{
					flag = false;
					break;
				}
			}
			if (flag == true)
			{
				temp[p++] = tmp_ind;
			}
		}
		diag_cd = p;
		//qsort(&temp[0], m, sizeof(INTTYPE), compare);
	//	for (int i = 0; i < p; i++)
	//	{
	//		cout << temp[i] << endl;;
		//}
		delete[] temp;
	}
//rows
	void CooMatrix::SortWithVectorsRow()
	{
		INTTYPE k = 0;
		vector < vector<INTTYPE> > vec_col_ind;
		vec_col_ind.resize(N);
		vector< vector<FPTYPE> > vec_val;
		vec_val.resize(N);
		for (INTTYPE i = 0; i < NNZ; i++)
		{
			vec_col_ind[row_ind[i]].push_back(col_ind[i]);
			vec_val[row_ind[i]].push_back(val[i]);
		}
		for (INTTYPE i = 0; i < N; i++)
		{
			for (INTTYPE j = 0; j < vec_col_ind[i].size(); j++)
			{
				row_ind[k] = i;
				col_ind[k] = vec_col_ind[i][j];
				val[k] = vec_val[i][j];
				k++;
			}
		}
	}
//Cols sort
//Need if matrix is symmetric for column coo storage (as default for general matrix)
void CooMatrix::SortWithVectorsCol()
	{
		INTTYPE k = 0;
		vector < vector<INTTYPE> > vec_row_ind;
		vec_row_ind.resize(N);
		vector< vector<FPTYPE> > vec_val;
		vec_val.resize(N);

		for (INTTYPE i = 0; i < NNZ; i++)
		{
			vec_row_ind[col_ind[i]].push_back(row_ind[i]);
			vec_val[col_ind[i]].push_back(val[i]);
		}
		for (INTTYPE i = 0; i < N; i++)
		{
			for (INTTYPE j = 0; j < vec_row_ind[i].size(); j++)
			{
				col_ind[k] = i;
				row_ind[k] = vec_row_ind[i][j];
				val[k] = vec_val[i][j];
				k++;
			}
		}
	}
//Characteristics calculation
	void CooMatrix::calculate_params()
	{
		FPTYPE temp;
		qsort(nnz_row, N, sizeof(INTTYPE),compare);
		min_nnz_per_row = nnz_row[0];
		max_nnz_per_row = nnz_row[N-1];
		avg_nnz_per_row = NNZ/N;
		FPTYPE disp = (FPTYPE)NNZ/(FPTYPE)(N + 1);
		for (int i=0; i<N;i++)
		{
			sigma+=pow(((FPTYPE)nnz_row[i]-disp),2.0);
		}
		sigma/=(FPTYPE)N+1;
		temp = sqrt(sigma);
		sigma =temp;
	}