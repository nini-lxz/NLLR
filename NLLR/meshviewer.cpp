#include <iostream>
#include "MeshInterface.h"

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		printf("You must enter 5 arguments.\n");
		exit(-1);
	}
	else
	{
		string input_path = argv[1];
		double nsig = stod(argv[2]);
		int numSimilarPatch = stoi(argv[3]);
		int vertexIter = stoi(argv[4]);

		printf("input noisy model path is: %s\n", argv[1]);
		printf("sigma_m = %0.3f\n", nsig);
		printf("N_k = %d\n", numSimilarPatch);
		printf("v_iter = %d\n\n", vertexIter);

		MyMeshInit_Guided_NLLR_global_denoising(input_path, nsig, numSimilarPatch, vertexIter);
	}	
	return 0;
}