// TAUCS.h: interface for the TAUCS library functions.
//
//////////////////////////////////////////////////////////////////////

#ifndef _LINEAR_SOLVER_TAUCS_H_INCLUDED_
#define _LINEAR_SOLVER_TAUCS_H_INCLUDED_

// Factorization Options
#define TAUCS_FACTOR_OPTION_LLT						0	// Cholesky (Symmetric and Lower Storage)
#define TAUCS_FACTOR_OPTION_LU						1	// LU (Asymmetric)
#define TAUCS_FACTOR_OPTION_LDLT					2	// LDL^T
#define TAUCS_FACTOR_OPTION_MF						3	// Multi-frontal
#define TAUCS_FACTOR_OPTION_LL						4	// Left-looking
#define TAUCS_FACTOR_OPTION_LLT_REUSE_SYMBOLIC		5	// Cholesky - Reusing Symbolic Factorization

class CSpMatrix;

// TAUCS Computation
class CTAUCS
{
public:
	// Factorization
	// 1. mat: Input
	// 2. f: Output
	// 3. bSym: Whether matrix is symmetric or not	
	// 4. nOpt: Method
	static bool Factor(CSpMatrix &mat, void *&f, const bool bSym, const int nOpt);

	// Solving Ax = b with Factorization Result
	// 1. f: Factorization
	// 2. x: Variables
	// 3. b: Right Hand Side
	// 4. rhsNum: No. of Right Hand Side
	static bool Solve(void *f, double *x, double *b, const int rhsNum);

	// Solving Ax = b
	// 1. matA: Coefficient Matrix
	// 2. x: Variables
	// 3. b: Right Hand Side
	// 4. bSym: Whether matrix A is symmetric or not	
	// 5. rhsNum: No. of Right Hand Side
	// 6. nOpt: Method
	static bool Solve(CSpMatrix &matA, double *x, double *b, const bool bSym, const int rhsNum, const int nOpt);

	// De-allocation
	// 1. f: Factorization
	static void Delete(void *&f);

	// Allocation
	// 1. mat: Matrix to be Allocated or De-allocated
	// 2. rowNum, colNum: No. of Rows and Columns
	// 3. nnzNum: No. of Non-zero Entries
	// 4. bSym: Whether matrix is symmetric or not
	// 5. temp: Intermediate Matrix for Conversion
	// 6. f: Factorization
/*	static bool Create(taucs_ccs_matrix *&mat, const int rowNum, const int colNum, const int nnzNum, const bool bSym);
	static bool Create(taucs_ccs_matrix *&mat, const CMat4TAUCS &temp, const bool bSym);
	static void Delete(taucs_ccs_matrix *mat);
	static void Delete(void *&f);

	// Factorization
	// 1. matA: Input
	// 2. f: Output
	// 3. nOpt: Method
	static bool Factor(taucs_ccs_matrix *matA, void *&f, const taucs_factor_option nOpt);

	// Solving Ax = b with Factorization Result
	// 1. matA: Input
	// 2. f: Factorization
	// 3. matX: Variables
	// 4. matB: Right Hand Side
	// 5. rhsNum: No. of Right Hand Side
	static bool Solve(taucs_ccs_matrix *matA, void *f, void *matX, void *matB, const int rhsNum);

	// Solving Ax = b
	// 1. matA: Input
	// 2. matX: Variables
	// 3. matB: Right Hand Side
	// 4. rhsNum: No. of Right Hand Side
	// 5. nOpt: Factorization Method
	static bool Solve(taucs_ccs_matrix *matA, void *matX, void *matB, const int rhsNum, const taucs_factor_option nOpt);*/
};

#endif	// _LINEAR_SOLVER_TAUCS_H_INCLUDED_
