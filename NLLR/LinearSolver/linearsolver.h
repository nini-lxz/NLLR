#pragma once

//#define _SUPERLU

#include "SpMatrix.h"
#include "SuperLU.h"
#include "TAUCS.h"
#include <time.h>

#ifdef _DEBUG
#pragma comment(lib, "../LinearSolver/LinearSolverLib_D.lib")
#else
#pragma comment(lib, "../LinearSolver/LinearSolverLib_R.lib")
#endif

class LinearSolver
{
public:
	static bool FactorA(CSpMatrix &matA, void *&f)
	{
		bool e=true;
		clock_t time = clock();
		//printf("%d %d %d\n",matA.m_rowNum, matA.m_colNum, matA.m_nnzNum);
	#ifdef _SUPERLU
		e= CSuperLU::Factor(matA, f);
	#else
		//FILE *fp = fopen("A.txt", "w");
		//for (int col=0; col<matA.m_colNum; col++){
		//	for (int j = 0; j < matA.m_rowCols[col].m_num; j++){
		//		int row = matA.m_rowCols[col].m_ids[j];
		//		double value = matA.m_rowCols[col].m_vals[j];
		//		fprintf(fp, "%d %d %.10f\n", row+1, col+1, value);
		//	}
		//}
		//fclose(fp);

		e= CTAUCS::Factor(matA, f, true, TAUCS_FACTOR_OPTION_LLT);


	#endif
		//printf("Factorization time: %dms\n", clock()-time);
		return e;
	};

	static bool SolveX(void *f, double *x, double *b)
	{
		bool e =1;
		clock_t time = clock();
	#ifdef _SUPERLU
		e= CSuperLU::Solve(f, x, b, 1);
	#else
		//FILE *fp = fopen("b.txt", "w");
		//for (int i=0; i<tetraNum*3; i++)
		//	fprintf(fp, "%.10f\n", b[i]);
		//fclose(fp);
		
		e= CTAUCS::Solve(f, x, b, 1);
	#endif
		//printf("Solving time: %dms\n", clock()-time);
		return e;
	};

	static bool SolveX(CSpMatrix &matA, double *x, double *b)
	{
		bool e =1;
		clock_t time = clock();
	#ifdef _SUPERLU
		e= CSuperLU::Solve(matA, x, b, 1);
	#else
		//FILE *fp = fopen("b.txt", "w");
		//for (int i=0; i<tetraNum*3; i++)
		//	fprintf(fp, "%.10f\n", b[i]);
		//fclose(fp);
		
		e = CTAUCS::Solve(matA, x, b, true, 1, TAUCS_FACTOR_OPTION_LLT);
	#endif
		//printf("Solving time: %dms\n", clock()-time);
		return e;
	};


	static void DeleteF(void *f){
	#ifdef _SUPERLU
		CSuperLU::Delete(f);
	#else
		CTAUCS::Delete(f);
	#endif
	};
};