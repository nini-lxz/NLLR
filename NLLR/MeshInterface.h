#ifndef __WMQ_INTERFACE_H__
#define __WMQ_INTERFACE_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

///////////////////////////////////////////////////////////////////////////////
#include "mesh/extension/ExItems.h"
#include "mesh\extension/ExKernelT.cpp" 
#include "mesh\read_write/read_write.cpp" 
#include <string>
#include <ctime>

using namespace std;

typedef MeshN::ExKernelT<MeshN::ExItems>   MyMesh;
typedef MeshN::ReaderWriterT<MyMesh>       Reader;

typedef MeshN::ExKernelT<MeshN::ExItems>   GroundMesh;
typedef MeshN::ReaderWriterT<GroundMesh>       GroundReader;

GroundMesh  gMesh;
GroundReader gReader(&gMesh);
void MyMeshInit_Guided_NLLR_global_denoising(string input_path, double nsig, int numSimilarPatch, int vertexIter)
{
	//-------------read noisy model------------
	MyMesh  mesh;
	Reader  reader(&mesh);
	time_t start = 0, end = 0;
	reader.off_reader(input_path.c_str());
	mesh.meshInit();

	cout << "Start processing....." << endl;
	
	std::vector<MyMesh::Normal> updatedVertexNormals; //store filtered vertex normals of last iteration, xz
	std::vector<MyMesh::Normal> bilateralGuidanceNormals;
	//fixed parameters
	int candidateNeighbSize = 10; //8; 10
	int K = 50;

	//////////////////////////step 1: non-local similar vertex selection/////////////////////////////////////
	mesh.bilateralNormalFilter(bilateralGuidanceNormals, nsig);
	mesh.patchNormalAlignment(bilateralGuidanceNormals, K);
	start = clock();
	mesh.getRegionNormalCovarianceMat_new(bilateralGuidanceNormals, K);
	mesh.findSimilarPatch3(candidateNeighbSize, numSimilarPatch);
	end = clock();
	cout << "[Non-local Similar Vertex Selection] running time: " 
		<< double(end - start) / CLOCKS_PER_SEC 
		<< " s" << endl;

	/////////////////////////step 2: low rank recover////////////////////////////////////////////
	double c1 = 2.8*sqrt(2);
	start = clock();
	mesh.low_rank_recover(numSimilarPatch, c1, nsig, K);
	end = clock();
	cout << "[Low-rank Recovery] running time: " 
		<< double(end - start) / CLOCKS_PER_SEC 
		<< " s" << endl;

	///////////////////////////////step3: vertex updating////////////////////////////
	mesh.compute_normals_NLLR();
	start = clock();
	mesh.adjustVertexCoord(vertexIter);
	end = clock();
	cout << "[Vertex Updating] running time: " 
		<< double(end - start) / CLOCKS_PER_SEC << " s" << endl;
	mesh.meshInit();
	mesh.calcUpdatedVertexNormal(updatedVertexNormals);

	//save our intermediate result
	char filename[200];
	sprintf(filename, "result_sig=%0.3f_numSimilarPatch=%d_vertexIter=%d.off",
		    nsig, numSimilarPatch, vertexIter);
	mesh.output_to_file(filename);
	
	cout << "Ending of the process!!!" << endl;
}

#endif // 