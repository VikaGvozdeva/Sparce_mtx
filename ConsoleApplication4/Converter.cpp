#define _CRT_SECURE_NO_WARNINGS

#include "Converter.h"
void Converters::qs(INTTYPE* s_arr, INTTYPE* _s_arr, int first, int last)
{
	int i = first, j = last, x = s_arr[(first + last) / 2], tmp = 0;

	do {
		while (s_arr[i] > x) i++;
		while (s_arr[j] < x) j--;

		if (i <= j) {
			if (s_arr[i] < s_arr[j])
			{
				tmp = s_arr[i];
				s_arr[i] = s_arr[j];
				s_arr[j] = tmp;

				tmp = _s_arr[i];
				_s_arr[i] = _s_arr[j];
				_s_arr[j] = tmp;

			}
			i++;
			j--;
		}
	} while (i <= j);

	if (i < last)
		qs(s_arr, _s_arr, i, last);
	if (first < j)
		qs(s_arr, _s_arr, first, j);
}

//void SymCOOToCOO(const COOMatrix &Mtx, COOMatrix &Matrix);
	void Converters::COOToCRS(const CooMatrix &Mtx, CRSMatrix &Matrix)
	{
		Matrix.DeleteAndAllocateMemory(Mtx.NNZ, Mtx.N, Mtx.isSym);
		int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0, n = 0, m = 0;
		NNZ = Mtx.NNZ;
		N = Mtx.N;
		vector< vector<int> > col_ind_;
		col_ind_.resize(N);
		vector< vector<double> > val_;
		val_.resize(N);

		for (i = 0; i < NNZ; i++)
		{
			val_[Mtx.row_ind[i]].push_back(Mtx.val[i]);
			col_ind_[Mtx.row_ind[i]].push_back(Mtx.col_ind[i]);
			Matrix.row_ptr[Mtx.row_ind[i] + 1]++;
		}

		for (k = 2; k < N + 1; k++)
		{
			Matrix.row_ptr[k] += Matrix.row_ptr[k - 1];
		}

		for (m = 0; m < N; m++) {
			k = 0;
			for (n = Matrix.row_ptr[m]; n < Matrix.row_ptr[m + 1]; n++)
			{
				Matrix.val[n] = val_[m][k];
				Matrix.col_ind[n] = col_ind_[m][k];
				k++;
			}
		}
	}

	void Converters::COOToCCS(const CooMatrix &Mtx, CCSMatrix &Matrix)
	{
		Matrix.DeleteAndAllocateMemory(Mtx.NNZ, Mtx.N);

		int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0;
		NNZ = Mtx.NNZ;
		N = Mtx.N;

		for (i = 0; i < NNZ; i++)
		{
			Matrix.val[i] = Mtx.val[i];
			Matrix.row_ind[i] = Mtx.row_ind[i];
		}

		for (j = 0; j < NNZ; j++)
		{
			tmp_ind = Mtx.col_ind[j];
			Matrix.col_ptr[++tmp_ind]++;
		}
		while (j < NNZ);

		for (k = 2; k < N + 1; k++)
		{
			Matrix.col_ptr[k] += Matrix.col_ptr[k - 1];
		}
	}

	void Converters::COOToCD(const CooMatrix &Mtx, CDMatrix &Matrix)
	{

		Matrix.DeleteAndAllocateMemory(Mtx.NNZ, Mtx.N, Mtx.diag_cd);
		int i = 0, j = 0, n, p = 0, k = 0, l = 0, NNZ = 0, N = 0, diag_ind = 0, B = 0, tmp_ind = 0, m = 0, diag_numb = 0;
		NNZ = Matrix.NNZ;
		N = Matrix.N;
		diag_numb = Matrix.B;
//		memset(Matrix.val, 0, N * diag_numb * sizeof(FPTYPE));
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < diag_numb; j++)
			{
				Matrix.val[i][j] = 0;
			}
		}
		memset(Matrix.diag, 0, diag_numb * sizeof(INTTYPE));
		Matrix.FillDiagArray(Mtx);
		vector< vector<FPTYPE> > val_;
		val_.resize(diag_numb);

		vector< vector<INTTYPE> > row_ind_;
		row_ind_.resize(diag_numb);

		for (int i = 0; i < diag_numb; i++)
		{
			val_[i].push_back(Matrix.diag[i]);
		}

		for (i = 0; i < NNZ; i++)
		{
			tmp_ind = Mtx.col_ind[i] - Mtx.row_ind[i];
			for (j = 0; j < diag_numb; j++)
			{
				if (tmp_ind == val_[j][0])
				{
					val_[j].push_back(Mtx.val[i]);
					row_ind_[j].push_back(Mtx.row_ind[i]);
			}
		}
		INTTYPE row_ind;
		for (m = 0; m < diag_numb; m++) {
			for (n = 1; n < val_[m].size(); n++)
			{
				row_ind = row_ind_[m][n - 1];
				Matrix.val[row_ind][m] = val_[m][n];
				}
			}
		}
		val_.clear();
		row_ind_.clear();
	}

	void Converters::COOToJD(const CooMatrix &Mtx, JDMatrix &Matrix)
	{
		Matrix.DeleteAndAllocateMemory(Mtx.NNZ, Mtx.N, Mtx.maxval_jd);
		int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0, MaxNNZ = 0;
		NNZ = Matrix.NNZ;
		N = Matrix.N;
		MaxNNZ = Matrix.MaxNNZ;
		vector< vector<INTTYPE> > col_ind_;
		col_ind_.resize(N);
		vector< vector<FPTYPE> > val_;
		val_.resize(N);
		INTTYPE *nnz_row = new INTTYPE[N];
		INTTYPE *nnz_col = new INTTYPE[N];
		INTTYPE *nnz_row_ind = new INTTYPE[N];

		for (i = 0; i < N; i++)
		{
			nnz_row[i] = Mtx.nnz_row[i];
			nnz_row_ind[i] = i;
			nnz_col[i] = 0;
		}
		//Converters::qs(nnz_row, nnz_row_ind, 0, N - 1);
		Converters::qs(nnz_row, nnz_row_ind, 0, N );

		for (i = 0; i < N; i++)
		{
			Matrix.perm[i] = nnz_row_ind[i];
		}
		for (i = 0; i < NNZ; i++)
		{
			val_[Mtx.row_ind[i]].push_back(Mtx.val[i]);
			col_ind_[Mtx.row_ind[i]].push_back(Mtx.col_ind[i]);
		}
		int sum = 0;
		for (int i = 0; i < N; i++)
		{
			int temp_nnz_col = val_[i].size();
			sum += temp_nnz_col;
			//int temp_nnz_col = nnz_row[Matrix.perm[i]];
			for (int k = 0; k < temp_nnz_col; k++)
			{
				nnz_col[k]++;
			}
		}
		int a = sum;
		sum = 0;
		int l = 0;
		int m = 0;
		for (i=0; i<MaxNNZ;i++)
		{
			//for (j = 0; j < nnz_row[Matrix.perm[i]];j++)
			sum += nnz_col[i];
			for (j = 0; j < nnz_col[i]; j++)
			{

				Matrix.jdiag[l] = val_[Matrix.perm[j]][i];
				Matrix.col_ind[l] = col_ind_[Matrix.perm[j]][i];
				l++;
			}
		}
		int b = sum;
		Matrix.jd_ptr[0] = 0;
		for (i = 1; i < MaxNNZ + 1; i++)
		{
			Matrix.jd_ptr[i] = nnz_col[i - 1];
		}

		for (i = 2; i < MaxNNZ + 1; i++)
		{
			Matrix.jd_ptr[i] += Matrix.jd_ptr[i - 1];
		}

		delete[] nnz_row;
		delete[] nnz_col;
		delete[] nnz_row_ind;

	}

 	void Converters::COOToSL(const CooMatrix &Mtx, SLMatrix &Matrix)
 	{
		Matrix.DeleteAndAllocateMemory(Mtx.NNZ, Mtx.N, Mtx.diag_elem, Mtx.isSym);
 		int j = 0, k = 0, l = 0, tmp_ind = 0, m = 0, p = 0, pos = 0;
 		int h = 0, o = 0;
 		int NNZ = Mtx.NNZ;
 		int N = Mtx.N;
 		vector< vector<FPTYPE> > vec;
 		vec.resize(N);

 		INTTYPE* elem_before_diag = new INTTYPE[N];
 		for (int i = 0; i < N; i++)
 		{
 			elem_before_diag[i] = 0;
 		}

		if (Mtx.isSym == 1)
		{
			for (int i = 0; i < NNZ; i++)
			{
				if (Mtx.col_ind[i] < Mtx.row_ind[i])
				{
					vec[Mtx.row_ind[i]].push_back(Mtx.val[i]);
					elem_before_diag[Mtx.row_ind[i]]++;
					Matrix.jptr[pos] = Mtx.col_ind[i];
					pos++;
				}
				else if (Mtx.col_ind[i] == Mtx.row_ind[i])
				{
					Matrix.adiag[Mtx.col_ind[i]] = Mtx.val[i];
				}
			}
		}
 		else if (Mtx.isSym == 0)
 		{
	 		for (int i = 0; i < NNZ; i++)
	 			{
	 				if (Mtx.col_ind[i] < Mtx.row_ind[i])
	 				{
	 					vec[Mtx.row_ind[i]].push_back(Mtx.val[i]);
	 					elem_before_diag[Mtx.row_ind[i]]++;
	 				}
	 				else if (Mtx.col_ind[i] > Mtx.row_ind[i])
	 				{
	 					Matrix.autr[pos] = Mtx.val[i];
	 					Matrix.jptr[pos] = Mtx.row_ind[i];
	 					pos++;
	 				}
	 				else if (Mtx.col_ind[i] == Mtx.row_ind[i])
	 				{
	 					Matrix.adiag[Mtx.col_ind[i]] = Mtx.val[i];
	 				}			
	 			}

 		}
 		Matrix.iptr[0] = 0;
 		for (int i = 0; i < N; i++)
 		{
 			if (elem_before_diag[i] != 0)
 			{
 				Matrix.iptr[i+1] = elem_before_diag[i];
 				for (int j = 0; j < elem_before_diag[i]; j++)
 				{
 					Matrix.altr[m] = vec[i][j];
 					m++;
 				}
 			}
 		}

 		Matrix.iptr[1] += Matrix.iptr[2];
 		for (p = 2; p < N+1; p++)
 		{
 			Matrix.iptr[p] += Matrix.iptr[p - 1];
 		}
 	
 		delete[] elem_before_diag;
 }
	