#define _CRT_SECURE_NO_WARNINGS

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <mkl.h>
#include <omp.h>
#include "Converter.h"
#include "CRS.h"
#include "COO.h"
#include "CCS.h"
#include "JD.h"
#include "CD.h"
#include "SL.h"

using namespace std;
#define TSC 2200000000.0
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
#else
typedef float FPTYPE;
#define ZERO 0.0f
#endif


FPTYPE SearchMax_double(FPTYPE* arr, INTTYPE N)
{
	int i;
	FPTYPE max_arr = arr[0];
	for (i = 1; i < N; i++)
	{
		if (max_arr < arr[i])
			max_arr = arr[i];
	}
	return max_arr;
}
FPTYPE CheckCorrectness(FPTYPE* my_mult, FPTYPE* mkl_mult, INTTYPE N)
{
	INTTYPE i;
	FPTYPE res;
	FPTYPE* arr_abs = new FPTYPE[N];
	for (i = 0; i < N; i++)
		arr_abs[i] = abs(my_mult[i] - mkl_mult[i]);
	res = SearchMax_double(arr_abs, N);
	delete[] arr_abs;
	return res;
}


/*FPTYPE CheckError(INTTYPE N, char* fileMKL, char* fileMult)
{
	FPTYPE res;
	FPTYPE *my_mult = new FPTYPE[N];
	FPTYPE *mkl_mult = new FPTYPE[N];
	FPTYPE* arr_abs = new FPTYPE[N];
	FILE *MKL = fopen(fileMKL, "rb");
	if (MKL == NULL)
	{
		printf("Error opening file");
	}
	fread(mkl_mult, sizeof(FPTYPE), N, MKL);

	FILE *MULT = fopen(fileMult, "rb");
	if (MULT == NULL)
	{
		printf("Error opening file");
	}
	fread(my_mult, sizeof(FPTYPE), N, MULT);

	for (int i = 0; i < N; i++)
		arr_abs[i] = abs(my_mult[i] - mkl_mult[i]);
	res = SearchMax_double(arr_abs, N);
	delete[] arr_abs;
	delete[] mkl_mult;
	delete[] my_mult;
	return res;

}
*/

char* getFileNameFromPath(char* path)
{
	for (size_t i = strlen(path) - 1; i >= 0; i--)
	{
		if ((path[i] == '/') || (path[i] == '\\'))
		{
			return &path[i + 1];
		}
	}
	printf("%s \n", path);
	return path;
}

void writeResultVector(FPTYPE* vec, INTTYPE N, char* filename)
{
	FILE *result = NULL;
	result = fopen(filename, "wb");
	if (result == NULL)
	{
		printf("Error opening file");
	}
	fwrite(vec, sizeof(FPTYPE), N, result);
	fclose(result);
}

long long ReadTSC() { // Returns time stamp counter
	int dummy[4]; // For unused returns
	volatile int DontSkip; // Volatile to prevent optimizing
	long long clock; // Time
	__cpuid(dummy, 0); // Serialize
	DontSkip = dummy[0]; // Prevent optimizing away cpuid
	clock = __rdtsc(); // Read time
	return clock;
}
int main(int argc, char** argv)
{

	char* fileName;
	char *act;
	char *fileInput;
	char *fileOutput;
	char *typeMatrix;
	char *fileInputCoo;
	char *fileMKLRes;
	char *fileMultRes;
	char *fn_pnt;
	char *threads;

	FPTYPE* result_crs;
	FPTYPE* result_mkl;
	FPTYPE* result_ccs;
	FPTYPE* result_jd;
	FPTYPE* result_cd;
	FPTYPE* result_sl;

	FILE *fp = stdout;

	CooMatrix Matrix(0, 0, 0);
	CRSMatrix CRS(0, 0, 0);
	CCSMatrix CCS(0, 0);
	JDMatrix JD(0, 0, 0);
	CDMatrix CD(0, 0, 0);
	SLMatrix SL(0, 0, 0, 0);

	fileName = new char[300];
	fileInput = new char[300];
	fileMKLRes = new char[300];
	fileMultRes = new char[300];
	fileInputCoo = new char[300];
	fileOutput= new char[300];
	typeMatrix = new char[300];
	act = new char[3];
	threads = new char[2];

	long long startTime, endTime;

	strcpy(fileName, argv[1]);
	fprintf(fp, "Matrix file name: %s\n\n", fileName);
	char* name = new char[strlen(fileName) - 4];

	fn_pnt = getFileNameFromPath(fileName);
	strncpy(name, fn_pnt, strlen(fn_pnt) - 4);

	name[strlen(fn_pnt) - 4] = '\0';
	strcpy(act, argv[2]);
	strcpy(typeMatrix, argv[3]);
//	strcpy(threads, argv[4]);
	//int thr = atoi(threads);
	//int thr = 1;
	//cout << thr;
	//omp_set_num_threads(16);
	omp_set_dynamic(0); 
	omp_set_num_threads(1);
	if (strcmp("COO", typeMatrix) == 0)
	{
		//int a = omp_get_num_threads();
		//cout << "count of threads -" << omp_get_num_threads() << endl;
		sprintf(fileInput, "COO_%s.bin", name);
		sprintf(fileOutput, "COO_%s_res.dt", name);

		if (strcmp("wr", act) == 0)
		{
			fp = fopen(fileOutput, "a");
			cout <<"before read"<<endl;
			startTime =  ReadTSC();
			Matrix.ReadMatrix(fileName);
			endTime =  ReadTSC();
			cout <<"after read"<<endl;

			Matrix.SortWithVectorsCol();
			//cout << Matrix << endl;
			fprintf(fp, "Time read matrix:  \t\t%lf\n", (double)(startTime - endTime) / TSC);
			Matrix.WriteInBinaryFile(fileInput);
		/*	for (int i = 0; i < Matrix.N; i++)
				cout << Matrix.nnz_row[i] << " " << endl;*/
			fprintf(fp, "Matrix size: N %d NNZ %d Sym %d\n\n", Matrix.N, Matrix.NNZ, Matrix.isSym);
		}
		if (strcmp("r", act) == 0)
		{
			Matrix.ReadFromBinaryFile(fileInput);
			//for (int i = 0; i < Matrix.N; i++)
			//	cout << Matrix.nnz_row[i] << " " << endl;
			//cout << "NNZ- " << Matrix.NNZ << " N-" << Matrix.N << "Sym - " << Matrix.isSym << endl;
		}
		Matrix.calculate_params();
		fprintf(fp, "Matrix size: sigma -  %lf avg_nnz - %d min_nnz %d  max_nnz %d\n\n", Matrix.sigma, Matrix.avg_nnz_per_row, Matrix.min_nnz_per_row, Matrix.max_nnz_per_row);
		//cout << "NNZ- " << Matrix.NNZ << " N-" << Matrix.N << "Sym - " << Matrix.isSym << endl;
	//	cout << Matrix;
		//cout << "sigma- " << Matrix.sigma << " avg-" << Matrix.avg_nnz_per_row << "min - " << Matrix.min_nnz_per_row << " max-" << Matrix.max_nnz_per_row<< endl;
	}

	if (strcmp("CRS", typeMatrix) == 0)
	{	//omp_set_num_threads(16);

		sprintf(fileInput, "CRS_%s.bin", name);
		sprintf(fileInputCoo, "COO_%s.bin", name);
		sprintf(fileOutput, "CRS_%s_res.dt", name);
		fp = fopen(fileOutput, "a");
				cout << "count of threads -" << omp_get_num_threads() << endl;

		Matrix.ReadFromBinaryFile(fileInputCoo);
	//	cout << "NNZ- " << Matrix.NNZ << " N-" << Matrix.N << "Sym - "<<Matrix.isSym<<endl;
		fprintf(fp, "Matrix size: N %d NNZ %d\n\n", Matrix.N, Matrix.NNZ);
		FPTYPE* v = new FPTYPE[Matrix.N];
		memset(v, 1, sizeof(FPTYPE) * Matrix.N);

		if (strcmp("wr", act) == 0)
		{
			startTime = ReadTSC();
			Converters::COOToCRS(Matrix, CRS);
			endTime = ReadTSC();
			//cout << "TSC - " << endTime - startTime << endl;
			//fprintf(fp, "Time Convert in \n\nCRS:  \t\t%lf\n", (double)(startTime - endTime) / TSC);
			CRS.WriteInBinaryFile(fileInput);
		}
		if (strcmp("r", act) == 0)
		{
			CRS.ReadFromBinaryFile(fileInput);
		}

		result_crs = new FPTYPE[Matrix.N];

	//	if (Matrix.isSym == 0)
		//{
			startTime = ReadTSC();
			CRS.MatrixVectorMultCRS(v, CRS.N, result_crs);
			endTime = ReadTSC();
			fprintf(fp, "Time Matrix-Vector multiplication in \n\nCRS: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			printf("Time Matrix-Vector multiplication in \n\ntime: \t\t%lf\n", (double)(endTime - startTime) / TSC);

	/*	}
		else
		{
			cout << "sym" << endl;
			startTime = ReadTSC();
			CRS.SymMatrixVectorMultCRS(v, Matrix.N, result_crs);
			endTime = ReadTSC();
			fprintf(fp, "Time Matrix-Vector multiplication in \n\nCRS: \t\t%lf\n", (double)(endTime - startTime) / TSC);

		}*/
		sprintf(fileMultRes, "CRS_res_%s.bin", name);
		//sprintf(fileMKLRes, "MKL_res_%s.bin", name);
		//writeResultVector(result_crs, CRS.N, fileMultRes);
		//fprintf(fp, "Error \n\nMKL: \t\t%lf\n", CheckError(CRS.N, fileMKLRes, fileMultRes));*/
		FILE* res = fopen(fileMultRes, "wb");
			if (res == NULL)
			{
				printf("Error opening file");
			}
			else
			{
				fwrite(result_crs, sizeof(FPTYPE), CRS.N, res);
			}
	
			fclose(res);

		delete[] v;
		delete[] result_crs;
	}

	if (strcmp("CCS", typeMatrix) == 0)
{
		sprintf(fileInput, "CCS_%s.bin", name);
		sprintf(fileInputCoo, "COO_%s.bin", name);
		sprintf(fileOutput, "CCS_%s_res.dt", name);
		fp = fopen(fileOutput, "a");

		Matrix.ReadFromBinaryFile(fileInputCoo);
		//cout << "NNZ- " << Matrix.NNZ << " N-" << Matrix.N << endl;
		fprintf(fp, "Matrix size: N %d NNZ %d\n\n", Matrix.N, Matrix.NNZ);
		FPTYPE* v = new FPTYPE[Matrix.N];
		memset(v, 1, sizeof(FPTYPE) * Matrix.N);

		if (strcmp("wr", act) == 0)
		{
			startTime = ReadTSC();
			Converters::COOToCCS(Matrix, CCS);
			endTime = ReadTSC();
		//	cout << "TSC - " << endTime - startTime << endl;
		//	fprintf(fp, "Time Convert in \n\nCCS:  \t\t%lf\n", (double)(startTime - endTime) / TSC);
			CCS.WriteInBinaryFile(fileInput);
		}
		if (strcmp("r", act) == 0)
		{
			CCS.ReadFromBinaryFile(fileInput);
		}

		result_ccs = new FPTYPE[Matrix.N];
		result_mkl = new FPTYPE[Matrix.N];

		/*if (Matrix.isSym == 0)
		{*/	//omp_set_num_threads(1);

			startTime = ReadTSC();
			CCS.MatrixVectorMultCCS(v, CCS.N, result_ccs);
			endTime = ReadTSC();
			fprintf(fp, "Time Matrix-Vector multiplication in \n\nCCS: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			printf("Time Matrix-Vector multiplication in \n\ntime: \t\t%lf\n", (double)(endTime - startTime) / TSC);
	/*	}
		else
		{
			startTime = ReadTSC();
			CCS.SymMatrixVectorMultCCS(v, Matrix.N, result_ccs);
			endTime = ReadTSC();
			fprintf(fp, "Time Matrix-Vector multiplication in \n\nCCS: \t\t%lf\n", (double)(endTime - startTime) / TSC);
		}*/

		sprintf(fileMultRes, "CCS_res_%s.bin", name);
		sprintf(fileMKLRes, "MKL_res_%s.bin", name);
		//writeResultVector(result_ccs, CCS.N, fileMultRes);
		FILE* res = fopen(fileMKLRes, "rb");
			if (res == NULL)
			{
				printf("Error opening file");
			}
			else
			{
				fread(result_mkl, sizeof(FPTYPE), CCS.N, res);
			}
	
			fclose(res);
			fprintf(fp, "Error \n\nMKL: \t\t%lf\n", CheckCorrectness(result_ccs, result_mkl,CCS.N));
		
		delete[] result_mkl;
		delete[] result_ccs;
		delete[] v;
	}

	if (strcmp("MKL", typeMatrix) == 0)
	{
		sprintf(fileOutput, "MKL_%s_res.dt", name);
		fp = fopen(fileOutput, "a");
		sprintf(fileInput, "CRS_%s.bin", name);
		CRS.ReadFromBinaryFile(fileInput);
		fprintf(fp, "Matrix size: N %d NNZ %d\n\n", CRS.N, CRS.NNZ);
		FPTYPE* v = new FPTYPE[CRS.N];
		memset(v, 1, sizeof(FPTYPE) * CRS.N);
		for (int i = 0; i < CRS.N + 1; i++)
		{
			CRS.row_ptr[i]++;
		}
		for (int i = 0; i < CRS.NNZ; i++)
		{
			CRS.col_ind[i]++;
		}
		result_mkl = new FPTYPE[CRS.N];
		result_crs = new FPTYPE[CRS.N];
	//	if (CRS.isSym == 0)
		//{
			startTime = ReadTSC();
			mkl_dcsrgemv("N", &(CRS.N), CRS.val, CRS.row_ptr, CRS.col_ind, v, result_mkl);
			endTime = ReadTSC();
			fprintf(fp, "Time Matrix-Vector multiplication in \n\nMKL: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			printf("Time Matrix-Vector multiplication in \n\ntime: \t\t%lf\n", (double)(endTime - startTime) / TSC);

		//}
		//else
		//{
			/*startTime = ReadTSC();
			mkl_dcsrsymv("L", &(CRS.N), CRS.val, CRS.row_ptr, CRS.col_ind, v, result_mkl);
			endTime = ReadTSC();
			fprintf(fp, "Time Matrix-Vector multiplication in \n\nMKL: \t\t%lf\n", (double)(endTime - startTime) / TSC);*/
	//	}
		sprintf(fileMKLRes, "MKL_res_%s.bin", name);
		sprintf(fileMultRes, "CRS_res_%s.bin", name);
		//writeResultVector(result_mkl, CRS.N, fileMKLRes);
		FILE* res_mkl = fopen(fileMKLRes, "wb");
		if (res_mkl == NULL)
		{
			printf("Error opening file");
		}
		else
		{
			fwrite(result_mkl, sizeof(FPTYPE), CRS.N, res_mkl);
		}

		fclose(res_mkl);
		FILE* res = fopen(fileMultRes, "rb");
			if (res == NULL)
			{
				printf("Error opening file");
			}
			else
			{
				fread(result_crs, sizeof(FPTYPE), CRS.N, res);
			}
	
			fclose(res);

			fprintf(fp, "Error \n\nMKL: \t\t%lf\n", CheckCorrectness(result_crs, result_mkl,CRS.N));
		delete[] result_mkl;
		delete[] result_crs;
		delete[] v;

	}

	if (strcmp("JD", typeMatrix) == 0)
	{
		sprintf(fileInput, "COO_%s.bin", name);
		Matrix.ReadFromBinaryFile(fileInput);
		//Matrix.calculate_params();
	
		sprintf(fileOutput, "JD_%s_res.dt", name);
		fp = fopen(fileOutput, "a");
		fprintf(fp, "Matrix size: N %d NNZ %d\n\n", Matrix.N, Matrix.NNZ);
		//maxvalJD = Matrix.maxvalJDMatrix(Matrix);
		FPTYPE* v = new FPTYPE[Matrix.N];
		memset(v, 1, sizeof(FPTYPE) * Matrix.N);
		if (strcmp("wr", act) == 0)
		{
			//maxvalJD = Matrix.maxvalJDMatrix(Matrix);
			//cout << "first:" << endl;
			//for (int i = 0; i < Matrix.N; i++)
			//	cout << Matrix.nnz_row[i] << " " << endl;
			Matrix.maxvalJDMatrix(Matrix);
			////Matrix.calculate_params();
			//cout << "second:" << endl;
			//for (int i = 0; i < Matrix.N; i++)
			//	cout << Matrix.nnz_row[i] << " " << endl;
			startTime = ReadTSC();
			Converters::COOToJD(Matrix, JD);
			endTime = ReadTSC();
			//fprintf(fp, "Time Convert in \n\nJD: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			sprintf(fileInput, "JD_%s.bin", name);
			JD.WriteInBinaryFile(fileInput);
		}
		if (strcmp("r", act) == 0)
		{
			sprintf(fileInput, "JD_%s.bin", name);
			JD.ReadFromBinaryFile(fileInput);
		}

		result_jd = new FPTYPE[Matrix.N];
		result_mkl = new FPTYPE[Matrix.N];
	/*	if (Matrix.isSym == 0)
		{*/
			startTime = ReadTSC();
			JD.MatrixVectorMultJD(v, JD.N, result_jd);
			endTime = ReadTSC();
			fprintf(fp, "\nTime Matrix-Vector multiplication in \n\nJD: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			printf("Time Matrix-Vector multiplication in \n\ntime: \t\t%lf\n", (double)(endTime - startTime) / TSC);

		//}
		//else
		//{
		//	startTime = ReadTSC();
		//	JD.SymMatrixVectorMultJD(v, Matrix.N, result_jd);
		//	endTime = ReadTSC();
		//	fprintf(fp, "\nTime Matrix-Vector multiplication in \n\nJD: \t\t%lf\n", (double)(endTime - startTime) / TSC);
		//}

		sprintf(fileMultRes, "JD_res_%s.bin", name);
		sprintf(fileMKLRes, "MKL_res_%s.bin", name);
		FILE* res = fopen(fileMKLRes, "rb");
			if (res == NULL)
			{
				printf("Error opening file");
			}
			else
			{
				fread(result_mkl, sizeof(FPTYPE), JD.N, res);
			}
	
			fclose(res);
			fprintf(fp, "Error \n\nMKL: \t\t%lf\n", CheckCorrectness(result_jd, result_mkl,JD.N));
		
		delete[] v;
		delete[] result_jd;
		delete[] result_mkl;
	}

	if (strcmp("SL", typeMatrix) == 0)
	{

		sprintf(fileInput, "COO_%s.bin", name);
		Matrix.ReadFromBinaryFile(fileInput);
		sprintf(fileOutput, "SL_%s_res.dt", name);
		fp = fopen(fileOutput, "a");
		fprintf(fp, "Matrix size: N %d NNZ %d\n\n", Matrix.N, Matrix.NNZ);
		//SLdiag = Matrix.diagSLMatrix(Matrix);
		FPTYPE* v = new FPTYPE[Matrix.N];
		memset(v, 1, sizeof(FPTYPE) * Matrix.N);
		if (strcmp("wr", act) == 0)
		{
			Matrix.diagSLMatrix(Matrix);
			startTime = ReadTSC();
			Converters::COOToSL(Matrix, SL);
			endTime = ReadTSC();
			//fprintf(fp, "\nTime convert SL: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			sprintf(fileInput, "SL_%s.bin", name);
			SL.WriteInBinaryFile(fileInput);
		}
		if (strcmp("r", act) == 0)
		{
			sprintf(fileInput, "SL_%s.bin", name);
			SL.ReadFromBinaryFile(fileInput);
		}

		result_sl = new FPTYPE[Matrix.N];
		result_mkl = new FPTYPE[Matrix.N];

		startTime = ReadTSC();
		SL.MatrixVectorMultSL(v, SL.N, result_sl);
		endTime = ReadTSC();
		fprintf(fp, "Time Matrix-Vector multiplication in \n\nSL: \t\t%lf\n", (double)(endTime - startTime) / TSC);
		printf("Time Matrix-Vector multiplication in \n\ntime: \t\t%lf\n", (double)(endTime - startTime) / TSC);

		sprintf(fileMultRes, "SL_res_%s.bin", name);
		sprintf(fileMKLRes, "MKL_res_%s.bin", name);
		//writeResultVector(result_sl, SL.N, fileMultRes);
		FILE* res = fopen(fileMKLRes, "rb");
			if (res == NULL)
			{
				printf("Error opening file");
			}
			else
			{
				fread(result_mkl, sizeof(FPTYPE),SL.N, res);
			}
	
			fclose(res);
			fprintf(fp, "Error \n\nMKL: \t\t%lf\n", CheckCorrectness(result_sl, result_mkl,SL.N));
		
		delete[] v;
		delete[] result_sl;
		delete[] result_mkl;

	}
	if (strcmp("CD", typeMatrix) == 0)
	{

		sprintf(fileInput, "CD_%s.bin", name);
		sprintf(fileInputCoo, "COO_%s.bin", name);
		sprintf(fileOutput, "CD_%s_res.dt", name);
		fp = fopen(fileOutput, "a");
		Matrix.ReadFromBinaryFile(fileInputCoo);
		//CDdiag = Matrix.DiagCDMatrix(Matrix);
		Matrix.DiagCDMatrix(Matrix);
		//fp = fopen(fileOutput, "w");
		//cout << "cd_diag" << Matrix.diag_cd << endl;
		fprintf(fp, "Matrix size: N %d NNZ %d\n\n", Matrix.N, Matrix.NNZ);
		FPTYPE* v = new FPTYPE[Matrix.N];
		memset(v, 1, sizeof(FPTYPE) * Matrix.N);


		if (strcmp("wr", act) == 0)
		{
		//	cout << "before conv" << endl;
			//CD.FillDiagArray(Matrix);
			startTime = ReadTSC();
			Converters::COOToCD(Matrix, CD);
			endTime = ReadTSC();
		//	cout << "after conv" << endl;
	//		fprintf(fp, "Time Convert in \n\nCD: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			//sprintf(fileInput, "CD_%s.bin", name);
			CD.WriteInBinaryFile(fileInput);
		}

		if (strcmp("r", act) == 0)
		{
			//sprintf(fileInput, "CD_%s.bin", name);
			CD.ReadFromBinaryFile(fileInput);
		}

		result_cd = new FPTYPE[Matrix.N];
		result_mkl = new FPTYPE[Matrix.N];


/*		if (Matrix.isSym == 0)
		{
		*/	startTime = ReadTSC();
			//cout << "before mult" << endl;
			CD.MatrixVectorMultCD(v, CD.N, result_cd);
			endTime = ReadTSC();
			//cout << "after mult" << endl;
			fprintf(fp, "\nTime Matrix-Vector multiplication in \n\nCD: \t\t%lf\n", (double)(endTime - startTime) / TSC);
			printf("Time Matrix-Vector multiplication in \n\ntime: \t\t%lf\n", (double)(endTime - startTime) / TSC);

		//}
		//else
		//{
		//	cout << "cd sym" << endl;
		//	startTime = ReadTSC();
		//	CD.SymMatrixVectorMultCD(v, Matrix.N, result_cd);
		//	endTime = ReadTSC();
		//	fprintf(fp, "\nTime Matrix-Vector multiplication in \n\nCD: \t\t%lf\n", (double)(endTime - startTime) / TSC);
		//}

		sprintf(fileMultRes, "CD_res_%s.bin", name);
		sprintf(fileMKLRes, "MKL_res_%s.bin", name);

		FILE* res = fopen(fileMKLRes, "rb");
			if (res == NULL)
			{
				printf("Error opening file");
			}
			else
			{
				fread(result_mkl, sizeof(FPTYPE), CD.N, res);
			}
	
			fclose(res);
			fprintf(fp, "Error \n\nMKL: \t\t%lf\n", CheckCorrectness(result_cd, result_mkl, CD.N));
		//	fclose(fp);
		delete[] v;
		delete[] result_cd;
		delete[] result_mkl;
	}
	delete[] fileName;
	delete[] fileInput;
	delete[] fileMKLRes;
	delete[] fileInputCoo;
	delete[] fileOutput;
	delete[] typeMatrix;
	delete[] act;
	delete[] name;
	delete[] fileMultRes;


	return 0;
}

