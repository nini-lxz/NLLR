#include "..\LinearSolver\LinearSolver.h" //this is all you need to include for using
//default, TAUCS is using, if SuperLU is needed, go into LinerSolver.h to enable #define _SUPERLU

/////////////////////////////////////////////////////
//usage
/////////////////////////////////////////////////////
	time=clock();
	printf("Factorizating......");
	void *F=NULL;
	if (!LinearSolver::FactorA(K, F)){
		printf("Cannot factorize F!\n"); delete[]x; delete[]b; return false;
	}
	printf("(%dms)\n", clock()-time);

	time=clock();
	printf("Solving......");
	if (!LinearSolver::SolveX(F, x, b)){
		printf("Cannot solve!\n"); delete[]x; delete[]b; return false;
	}
	printf("(%dms)\n", clock()-time);

	//free memory
	LinearSolver::DeleteF(F);