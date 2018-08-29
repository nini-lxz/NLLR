// SuperLU.h: interface for the CSuperLU class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _LINEAR_SOLVER_SUPER_LU_H_INCLUDED_
#define _LINEAR_SOLVER_SUPER_LU_H_INCLUDED_

class CSpMatrix;

class CSuperLU	// SuperLU Computation
{
public:
	// Factorization
	// 1. mat: Input
	// 2. f: Output
	static bool Factor(CSpMatrix &mat, void *&f);

	// Solving Ax = b with Factorization Result
	// 1. f: Factorization
	// 2. x: Variables
	// 3. b: Right Hand Side
	// 4. rhsNum: No. of Right Hand Side
	static bool Solve(void *f, double *x, double *b, const int rhsNum);

	// Solving Ax = b
	// 1. matA: Input
	// 2. x: Variables
	// 3. b: Right Hand Side
	// 4. rhsNum: No. of Right Hand Side
	static bool Solve(CSpMatrix &matA, double *x, double *b, const int rhsNum);

	// De-allocation
	// 1. f: Factorization
	static void Delete(void *&f);
};

#endif	// _LINEAR_SOLVER_SUPER_LU_H_INCLUDED_
