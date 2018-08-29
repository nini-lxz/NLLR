
#include "ExKernelT.h"
#include <stdio.h>
#include <iostream>
#include "algorithm"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/IterativeLinearSolvers>

using namespace Eigen;
using namespace std;


namespace MeshN { 

	///////////////////////////////////////////////////////////////////////////////
	// Implementation of member functions of ExKernelT
	/////////////////////////////////////////////////////////////////////////////// 
	template<class ExItems>
	ExKernelT<ExItems>::ExKernelT(): KernelT<ExItems>() {

		//kdtree_ = NULL;
		//ps_     = NULL;
		isNormal_ = false;
		isArea_ = false;

	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	ExKernelT<ExItems>::~ExKernelT(){

		//if(kdtree_ != NULL) delete kdtree_;

		//if(ps_ != NULL)     delete ps_;

	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::meshInit(){
		update_facet_normals();
		update_vertex_normals_max();
	}

	////////////////////////////////////////////////////////////////////////////////
	template <class ExItems>
	void ExKernelT<ExItems>::update_flag(void){
		FacetIterator fi = facet_begin();
		for (; fi != facet_end(); ++fi) {
			if (fi->status_.is_deleted()) continue;

			assert(fi->halfedge_handle_.is_valid());
			fi->flag_ = false;
		}

		VertexIterator vi = vertex_begin();
		for (; vi != vertex_end(); ++vi) {
			if (vi->status_.is_deleted()) continue;

			assert(vi->halfedge_handle_.is_valid());
			vi->flag_ = false;
		}
	}
	////////////////////////////get_original_vertex//////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::get_original_vertex(){
			VertexIterator vi(vertex_begin());
			for (; vi != vertex_end(); vi++){

				HalfedgeHandle& hh = vi->halfedge_handle_;
				VertexHandle&   vh = vertex_handle(hh);
				original_vertices_.push_back(vertex_ref(vh));
			}

		}

	///////////////////////get_original_facet//////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::get_original_facet(){
			FacetIterator fi(facet_begin());
			for (; fi < facet_end(); fi++){

				HalfedgeHandle& hh = fi->halfedge_handle_;
				FacetHandle&   fh = facet_handle(hh);
				original_facets_.push_back(facet_ref(fh));
			}

		}

	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::patchNormalAlignment(std::vector<Normal> guidanceNormals, int K)
	{
		VertexIterator vi = vertex_begin();
		for (vi; vi < vertex_end(); vi++)
		{
			vector<VertexHandle> neighborVertices;
			getNeighborVertices(vertex_handle(vi->halfedge_handle_), K, neighborVertices);
			//////////////////////////////////////////////////////////////////////////////////
			MatrixXd neighborVertexNormals;
			neighborVertexNormals.setZero(neighborVertices.size(), 3);
			double sumNx = 0; double sumNy = 0; double sumNz = 0;
			for (int i = 0; i < neighborVertices.size(); i++)
			{
				neighborVertexNormals(i, 0) = normal(neighborVertices[i]).data_[0];
				neighborVertexNormals(i, 1) = normal(neighborVertices[i]).data_[1];
				neighborVertexNormals(i, 2) = normal(neighborVertices[i]).data_[2];
				sumNx += neighborVertexNormals(i, 0);
				sumNy += neighborVertexNormals(i, 1);
				sumNz += neighborVertexNormals(i, 2);
			}
			double len = sqrt(sumNx*sumNx + sumNy*sumNy + sumNz*sumNz);
			double avergNx = sumNx/len;
			double avergNy = sumNy/len;
			double avergNz = sumNz/len;
			//cout << "average: " << avergNx << "  " << avergNy << "  " << avergNz << endl;
			Matrix3d coMatrix;
			coMatrix.setZero(3, 3);
			coMatrix = (neighborVertexNormals.transpose())*neighborVertexNormals;
			
			EigenSolver<Matrix3d> es(coMatrix);
			MatrixXd eigenVal = es.pseudoEigenvalueMatrix();
			//cout << "eigenVal= " << endl << eigenVal << endl;
			MatrixXd eigenVec = es.pseudoEigenvectors();
			//cout << "eigenVec= " << endl << eigenVec << endl;
			MatrixXd rotMat; rotMat.setZero(3, 3);
			float maxEigenVal = 0; float minEigenVal = 100000;
			int maxIdx = 0; int minIdx = 0;
			for (int i = 0; i < 3; i++)
			{
				float eigen_value = eigenVal(i, i);
				if (eigen_value >= maxEigenVal)
				{
					maxEigenVal = eigen_value;
					maxIdx = i;
				}
				if (eigen_value <= minEigenVal)
				{
					minEigenVal = eigen_value;
					minIdx = i;
				}
			}
			int middleIdx = 0;
			for (int i = 0; i < 3; i++)
			{
				if (i != maxIdx && i != minIdx)
				{
					middleIdx = i;
				}
			}

			///////////// get n //////////////
			Vector3f n;
			if ((avergNx*eigenVec(0, maxIdx)) > 0)
			{
				n << eigenVec(0, maxIdx), eigenVec(1, maxIdx), eigenVec(2, maxIdx);
			}
			else
			{
				n << eigenVec(0, maxIdx)*-1, eigenVec(1, maxIdx)*-1, eigenVec(2, maxIdx)*-1;
			}
			n.normalize();
			//cout << "n:" << endl << n << endl;
			//out << "v " << n(0)*2 << " " << n(1) * 2 << " " << n(2) * 2 << "\n";
			////////////// get b /////////////////
			Vector3f t;
			t << eigenVec(0, middleIdx), eigenVec(1, middleIdx), eigenVec(2, middleIdx);
			t.normalize();
			Vector3f b = n.cross(t);
			b.normalize();
			//cout << "b:" << endl << b << endl;
			//out << "v "<<b(0) * 2 << " " << b(1) * 2 << " " << b(2) * 2 << "\n";
			////////////// get t /////////////
			t = b.cross(n);
			t.normalize();
			////////////////// construct rot matrix /////////////
			rotMat(0, 0) = t(0);
			rotMat(0, 1) = t(1);
			rotMat(0, 2) = t(2);
			rotMat(1, 0) = b(0);
			rotMat(1, 1) = b(1);
			rotMat(1, 2) = b(2);
			rotMat(2, 0) = n(0);
			rotMat(2, 1) = n(1);
			rotMat(2, 2) = n(2);
			//cout << "rotMat= " << endl << rotMat << endl;
			(*vi).rotMat_ = rotMat;
			(*vi).patch_aligned_ = neighborVertexNormals*(rotMat.transpose());
			////////////////// guidance ///////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////
			MatrixXd neighborGuiVertexNormals;
			neighborGuiVertexNormals.setZero(neighborVertices.size(), 3);
			sumNx = 0; sumNy = 0; sumNz = 0;
			for (int i = 0; i < neighborVertices.size(); i++)
			{
				neighborGuiVertexNormals(i, 0) = guidanceNormals[neighborVertices[i].idx()].data_[0];
				neighborGuiVertexNormals(i, 1) = guidanceNormals[neighborVertices[i].idx()].data_[1];
				neighborGuiVertexNormals(i, 2) = guidanceNormals[neighborVertices[i].idx()].data_[2];
				sumNx += neighborGuiVertexNormals(i, 0);
				sumNy += neighborGuiVertexNormals(i, 1);
				sumNz += neighborGuiVertexNormals(i, 2);
			}
			len = sqrt(sumNx*sumNx + sumNy*sumNy + sumNz*sumNz);
			avergNx = sumNx / len;
			avergNy = sumNy / len;
			avergNz = sumNz / len;
			MatrixXd guiCoMatrix;
			guiCoMatrix.setZero(3, 3);
			guiCoMatrix = (neighborGuiVertexNormals.transpose())*neighborGuiVertexNormals;

			EigenSolver<Matrix3d> es_gui(guiCoMatrix);
			MatrixXd eigenVal_gui = es_gui.pseudoEigenvalueMatrix();
			//cout << "eigenVal= " << endl << eigenVal << endl;
			MatrixXd eigenVec_gui = es_gui.pseudoEigenvectors();
			//cout << "eigenVec= " << endl << eigenVec << endl;
			MatrixXd rotMat_gui; rotMat_gui.setZero(3, 3);
			float maxEigenVal_gui = 0; float minEigenVal_gui = 100000;
			int maxIdx_gui = 0; int minIdx_gui = 0;
			for (int i = 0; i < 3; i++)
			{
				float eigen_value = eigenVal_gui(i, i);
				if (eigen_value >= maxEigenVal_gui)
				{
					maxEigenVal_gui = eigen_value;
					maxIdx_gui = i;
				}
				if (eigen_value <= minEigenVal_gui)
				{
					minEigenVal_gui = eigen_value;
					minIdx_gui = i;
				}
			}
			int middleIdx_gui = 0;
			for (int i = 0; i < 3; i++)
			{
				if (i != maxIdx_gui && i != minIdx_gui)
				{
					middleIdx_gui = i;
				}
			}
			
			///////////// get n //////////////
			Vector3f n_gui;
			if ((avergNx*eigenVec_gui(0, maxIdx_gui)) > 0)
			{
				n_gui << eigenVec_gui(0, maxIdx_gui), eigenVec_gui(1, maxIdx_gui), eigenVec_gui(2, maxIdx_gui);
			}
			else
			{
				n_gui << eigenVec_gui(0, maxIdx_gui)*-1, eigenVec_gui(1, maxIdx_gui)*-1, eigenVec_gui(2, maxIdx_gui)*-1;
			}
			n_gui.normalize();

			////////////// get b /////////////////
			Vector3f t_gui;
			t_gui << eigenVec_gui(0, middleIdx_gui), eigenVec_gui(1, middleIdx_gui), eigenVec_gui(2, middleIdx_gui);
			t_gui.normalize();
			Vector3f b_gui = n_gui.cross(t_gui);
			b_gui.normalize();

			////////////// get t /////////////
			t_gui = b_gui.cross(n_gui);
			t_gui.normalize();

			////////////////// construct rot matrix /////////////
			rotMat_gui(0, 0) = t_gui(0);
			rotMat_gui(0, 1) = t_gui(1);
			rotMat_gui(0, 2) = t_gui(2);
			rotMat_gui(1, 0) = b_gui(0);
			rotMat_gui(1, 1) = b_gui(1);
			rotMat_gui(1, 2) = b_gui(2);
			rotMat_gui(2, 0) = n_gui(0);
			rotMat_gui(2, 1) = n_gui(1);
			rotMat_gui(2, 2) = n_gui(2);

			(*vi).patch_aligned_gui_ = neighborGuiVertexNormals*(rotMat_gui.transpose());
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getRegionNormalCovarianceMat_new(std::vector<Normal> guidanceNormals, int K)
	{
		VertexIterator vi = vertex_begin();
		for (vi; vi < vertex_end(); vi++)
		{
			MatrixXd patch_aligned = (*vi).patch_aligned_;
			double sumNx = 0; double sumNy = 0; double sumNz = 0;
			for (int i = 0; i < K; i++)
			{
				sumNx += patch_aligned(i, 0);
				sumNy += patch_aligned(i, 1);
				sumNz += patch_aligned(i, 2);
			}
			double len = sqrt(sumNx*sumNx + sumNy*sumNy + sumNz*sumNz);
			double avergNx = sumNx / len;
			double avergNy = sumNy / len;
			double avergNz = sumNz / len;
			for (int i = 0; i < K; i++)
			{
				patch_aligned(i, 0) -= avergNx;
				patch_aligned(i, 1) -= avergNy;
				patch_aligned(i, 2) -= avergNz;
			}
			MatrixXd coMatrix;
			coMatrix.setZero(3, 3);
			coMatrix = (patch_aligned.transpose())*patch_aligned;
			(*vi).coMatrix_ = coMatrix;

			////////////////////////guidance////////////////////////////////////////////////
			MatrixXd patch_aligned_gui = (*vi).patch_aligned_gui_;
			sumNx = 0; sumNy = 0; sumNz = 0;
			for (int i = 0; i < K; i++)
			{
				sumNx += patch_aligned_gui(i, 0);
				sumNy += patch_aligned_gui(i, 1);
				sumNz += patch_aligned_gui(i, 2);
			}
			len = sqrt(sumNx*sumNx + sumNy*sumNy + sumNz*sumNz);
			avergNx = sumNx / len;
			avergNy = sumNy / len;
			avergNz = sumNz / len;
			for (int i = 0; i < K; i++)
			{
				patch_aligned_gui(i, 0) -= avergNx;
				patch_aligned_gui(i, 1) -= avergNy;
				patch_aligned_gui(i, 2) -= avergNz;
			}
			MatrixXd guiCoMatrix;
			guiCoMatrix.setZero(3, 3);
			guiCoMatrix = (patch_aligned_gui.transpose())*patch_aligned_gui;
			(*vi).guiCoMatrix_ = guiCoMatrix;
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::findSimilarPatch3(int _candidateSize, int _t) 
	{
		VertexIterator vi = vertex_begin();
		for (; vi != vertex_end(); vi++)
		{
			MatrixXd refCoMatrix;
			refCoMatrix = (*vi).coMatrix_;
				
			MatrixXd refGuiCoMatrix;
			refGuiCoMatrix = (*vi).guiCoMatrix_;

			vector<VertexHandle> NeighborRing;
			getNeighborRing(vertex_handle(vi->halfedge_handle_), _candidateSize, NeighborRing);

			vector<vector<double> > distanceMatrix(NeighborRing.size(), vector<double>(2));
			for (int i = 0; i < NeighborRing.size(); i++)
			{
				MatrixXd candidateCoMatrix;
				candidateCoMatrix = vertex_ref(NeighborRing[i]).coMatrix_;

				MatrixXd candidateGuiCoMatrix;
				candidateGuiCoMatrix = vertex_ref(NeighborRing[i]).guiCoMatrix_;

				MatrixXd disMatrix = refCoMatrix - candidateCoMatrix;
				MatrixXd disMatrix2 = disMatrix.array()*disMatrix.array();
				double dis = disMatrix2.sum();			

				MatrixXd disGuiMatrix = refGuiCoMatrix - candidateGuiCoMatrix;
				MatrixXd disGuiMatrix2 = disGuiMatrix.array()*disGuiMatrix.array();
				double guiDis = disGuiMatrix2.sum();

				//---------------------------------------------------
				distanceMatrix[i][0] = dis*guiDis;
				distanceMatrix[i][1] = NeighborRing[i].idx();
			}

			sort(distanceMatrix.begin(), distanceMatrix.end());

			for (int i = 0; i < _t; i++)
			{
				(*vi).similarIdx_.push_back(distanceMatrix[i][1]);
			}
		}
	}

	template<class ExItems>
	double ExKernelT<ExItems>::computeSig(double nsig) {

		VertexIterator vi = vertex_begin();
		double meanDifFaceNormal = 0.0;
		FacetIterator fi = facet_begin();
		int i = 0;
		for (; fi != facet_end(); ++fi) {
			FacetHandle _fh = facet_handle(fi->halfedge_handle_);
			Normal normalDif = facet_ref(_fh).normal_ - facet_ref(_fh).oriNormal_;
			meanDifFaceNormal = meanDifFaceNormal + normalDif[0] * normalDif[0] + normalDif[1] * normalDif[1] + normalDif[2] * normalDif[2];
			i++;
			
		}
		meanDifFaceNormal = meanDifFaceNormal / i;
		cout << "meanFaceNormalDif=" << meanDifFaceNormal << endl;
		double sig = sqrt(abs(nsig*nsig - meanDifFaceNormal));
		return sig;
	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::updateNormalIterativeRegularization(double w) {
		FacetIterator fi = facet_begin();
		for (; fi != facet_end(); ++fi) {
			FacetHandle _fh = facet_handle(fi->halfedge_handle_);
			Normal nd = (facet_ref(_fh).oriNormal_ - facet_ref(_fh).normal_);
			nd[0] = w*nd[0];
			nd[1] = w*nd[1];
			nd[2] = w*nd[2];
			fi->normal_ = facet_ref(_fh).normal_ + nd;
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::low_rank_recover(int _t, double c1, double nsig, int K) 
	{
		int orderMethod = 3; // 1: random ordering; 2: ring-based random; 3: ring-based NPC
		
		VertexIterator vi = vertex_begin();
		for (; vi != vertex_end(); vi++)
		{
			// low-rank matrix recovery
			int maxRows = 0;
			HalfedgeHandle hh = vi->halfedge_handle_;
			VertexHandle vh = vertex_handle(vi->halfedge_handle_);
		
			maxRows = K;
			MatrixXd Bx, By, Bz; 
			Bx.setZero(maxRows, _t); //Bx: the x coordinate
			By.setZero(maxRows, _t); //By: the y coordinate
			Bz.setZero(maxRows, _t); //Bz: the z coordinate

			////////////// construct NPG matrix //////////////////
			vector<vector<VertexHandle> >vv(K, vector<VertexHandle>(_t)); //NPG VertexHandle matrix

			if (orderMethod == 1)
			{
				for (int i = 0; i < _t; i++)
				{
					vector<VertexHandle> vii;
					int dix = vertex_ref(vh).similarIdx_[i];
					getNeighborVertices(VertexHandle(dix), K, vii);
					for (int j = 0; j < K; j++)
					{
						vv[j][i] = vii[j];
					}
				}

				for (int i = 0; i < _t; i++)
				{
					int dix = vertex_ref(vh).similarIdx_[i];
					MatrixXd patch_aligned = vertex_ref(VertexHandle(dix)).patch_aligned_;
					for (int j = 0; j < K; j++)
					{
						Bx(j, i) = patch_aligned(j, 0);
						By(j, i) = patch_aligned(j, 1);
						Bz(j, i) = patch_aligned(j, 2);
					}
				}
				//// mess up the order totally
				srand(time(NULL));
				int index, tmp;
				VertexHandle tem_vh;
				for (int i = 0; i < _t; i++)
				{
					for (int j = 0; j < K; j++)
					{
						index = rand() % (K - j) + j;
						if (index != j)
						{
							tmp = Bx(j, i);
							Bx(j, i) = Bx(index, i);
							Bx(index, i) = tmp;

							tmp = By(j, i);
							By(j, i) = By(index, i);
							By(index, i) = tmp;

							tmp = Bz(j, i);
							Bz(j, i) = Bz(index, i);
							Bz(index, i) = tmp;

							tem_vh = vv[j][i];
							vv[j][i] = vv[index][i];
							vv[index][i] = tem_vh;
						}
					}
				}
			}
			
			if (orderMethod == 2)
			{
				vector<vector<int> >vv_ringIdx(K, vector<int>(_t));
				for (int i = 0; i < _t; i++)
				{
					vector<VertexHandle> vii;
					vector<int> vii_ringIdx;
					int dix = vertex_ref(vh).similarIdx_[i];
					getNeighborVertices_withRingIndex(VertexHandle(dix), K, vii, vii_ringIdx);
					for (int j = 0; j < K; j++)
					{
						vv[j][i] = vii[j];
						vv_ringIdx[j][i] = vii_ringIdx[j];
					}
					vii_ringIdx.clear();
				}

				for (int i = 0; i < _t; i++)
				{
					int dix = vertex_ref(vh).similarIdx_[i];
					MatrixXd patch_aligned = vertex_ref(VertexHandle(dix)).patch_aligned_;
					for (int j = 0; j < K; j++)
					{
						Bx(j, i) = patch_aligned(j, 0);
						By(j, i) = patch_aligned(j, 1);
						Bz(j, i) = patch_aligned(j, 2);
					}
				}

				// mess up the order within each ring
				srand(time(NULL));
				int index, tmp;
				int maxRingNum = 0;
				int oldVertexNum = 0;
				VertexHandle tem_vh;
				for (int i = 0; i < _t; i++)
				{
					// find the maximum ring number for each patch
					maxRingNum = vv_ringIdx[K - 1][i];
					oldVertexNum = 0;
					for (int j = 0; j <= maxRingNum; j++)
					{
						vector<VertexHandle> ringVertices;
						for (int k = 0; k < K; k++)
						{
							if (vv_ringIdx[k][i] == j)
								ringVertices.push_back(vv[k][i]);
						}
						//// mess up the order within the ring
						for (int k = 0; k < ringVertices.size(); k++)
						{
							index = rand() % (ringVertices.size() - k) + k;
							index = index + oldVertexNum;
							if (index != k + oldVertexNum)
							{
								tmp = Bx(k + oldVertexNum, i);
								Bx(k + oldVertexNum, i) = Bx(index, i);
								Bx(index, i) = tmp;

								tmp = By(k + oldVertexNum, i);
								By(k + oldVertexNum, i) = By(index, i);
								By(index, i) = tmp;

								tmp = Bz(k + oldVertexNum, i);
								Bz(k + oldVertexNum, i) = Bz(index, i);
								Bz(index, i) = tmp;

								tem_vh = vv[k + oldVertexNum][i];
								vv[k + oldVertexNum][i] = vv[index][i];
								vv[index][i] = tem_vh;
							}
						}
						oldVertexNum = oldVertexNum + ringVertices.size();
						ringVertices.clear();
					}
				}
			}
			
			if (orderMethod == 3)
			{
				for (int i = 0; i < _t; i++)
				{
					vector<VertexHandle> vii;
					int dix = vertex_ref(vh).similarIdx_[i];
					getNeighborVertices(VertexHandle(dix), K, vii);
					for (int j = 0; j < K; j++)
					{
						vv[j][i] = vii[j];
					}
				}

				for (int i = 0; i < _t; i++)
				{
					int dix = vertex_ref(vh).similarIdx_[i];
					MatrixXd patch_aligned = vertex_ref(VertexHandle(dix)).patch_aligned_;
					for (int j = 0; j < K; j++)
					{
						Bx(j, i) = patch_aligned(j, 0);
						By(j, i) = patch_aligned(j, 1);
						Bz(j, i) = patch_aligned(j, 2);
					}
				}
			}

			//-------------- low rank recovery --------------------
			MatrixXd Xx, Wx, Xy, Wy, Xz, Wz;
			// c++ implementation (version 1)
			//low_rank_ssc_C(Bx, c1, nsig, Xx, Wx);
			//low_rank_ssc_C(By, c1, nsig, Xy, Wy);
			//low_rank_ssc_C(Bz, c1, nsig, Xz, Wz);
			// c++ implementation (version 2)
			low_rank_ssc_C_TG(Bx, c1, nsig, Xx, Wx);
			low_rank_ssc_C_TG(By, c1, nsig, Xy, Wy);
			low_rank_ssc_C_TG(Bz, c1, nsig, Xz, Wz);

			//--------------- rotation back -----------------------
			MatrixXd Xx_recover, Xy_recover, Xz_recover;
			Xx_recover.setZero(K, _t); Xy_recover.setZero(K, _t); Xz_recover.setZero(K, _t);
			for (int i = 0; i < _t; i++)
			{
				int dix = vertex_ref(vh).similarIdx_[i];
				MatrixXd rotMat = vertex_ref(VertexHandle(dix)).rotMat_;
				MatrixXd rotMat_inverse = rotMat.transpose();  //transpose is equal to inverse
				MatrixXd mat;
				mat.setZero(K, 3);
				for (int j = 0; j < K; j++)
				{
					mat(j, 0) = Xx(j, i);
					mat(j, 1) = Xy(j, i);
					mat(j, 2) = Xz(j, i);
				}
				MatrixXd mat_align_back = mat*(rotMat_inverse.transpose());
				for (int j = 0; j < K; j++)
				{
					Xx_recover(j, i) = mat_align_back(j, 0);
					Xy_recover(j, i) = mat_align_back(j, 1);
					Xz_recover(j, i) = mat_align_back(j, 2);
				}
			}
			
			//----------- accumulated assign ------------------------
			for (int i = 0; i < _t; i++)
			{
				for (int j = 0; j < K; j++)
				{
					MathN::Vec3f newNor = MathN::Vec3f(Xx_recover(j, i), Xy_recover(j, i), Xz_recover(j, i));
					MathN::Vec3f newWei = MathN::Vec3f(Wx(j, i), Wy(j, i), Wz(j, i));
					vertex_ref(vv[j][i]).weights_ += newWei;
					vertex_ref(vv[j][i]).sumNormal_ += newNor;
				}
			}

		}
	}
		
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::low_rank_ssc_C(MatrixXd Y, double c1, double nsig, MatrixXd &output, MatrixXd &W) {
		int rows = Y.rows();
		int cols = Y.cols();
		MatrixXd mY;
		mY.setZero(rows, cols); // mean value (normal component) of those similar patches
		for (int i = 0; i < rows; i++)
		{
			double avg = 0;
			for (int j = 0; j < cols; j++)
				avg += Y(i, j);
			avg = avg / cols;

			for (int j = 0; j < cols; j++)
				mY(i, j) = avg;
		}
		Y = Y - mY;
		JacobiSVD<MatrixXd> svdX(Y, ComputeThinU | ComputeThinV);
		MatrixXd U0 = svdX.matrixU();
		MatrixXd Sigma = svdX.singularValues();
		MatrixXd V0 = svdX.matrixV();
		MatrixXd Sigma0 = Sigma.asDiagonal();
		//cout << Sigma << endl <<endl;
		
		int PatNum = cols;
		double eps = 0.00001;
		MatrixXd Temp = Sigma0.array() * Sigma0.array();

		Temp = Temp.array() - PatNum*nsig*nsig;
		
		for (int i = 0; i < Temp.rows(); i++)
		{
			for (int j = 0; j < Temp.cols(); j++)
			{
				if (Temp(i ,j) < 0)
				{
					Temp(i, j) = 0;
				}
			}

		}
		Temp = Temp.array().sqrt();
		MatrixXd thr = (c1*sqrt(PatNum)*(nsig*nsig)) / (Temp.array() + eps);
		Sigma0 = Sigma0.array() - thr.array();

		int r = 0;
		for (int i = 0; i < Sigma0.rows(); i++)
		{
			for (int j = 0; j < Sigma0.cols(); j++)
			{
				if (Sigma0(i, j) <= 0)
				{
					Sigma0(i, j) = 0;
				}
				else
				{
					r = r + 1;
				}
			}
		}
		if (r == 0)
		{
			output.setZero(rows, cols);
		}
		else
			output = U0.middleCols(0, r - 1)*Sigma0.middleCols(0, r - 1).middleRows(0, r - 1) *V0.middleCols(0, r - 1).transpose();
		double wei;
		//if (r == rows || r == 0)
		if (r == rows || r == 0)
		{
			wei = 1 / (double)rows;
		} 
		else
		{
			wei = (rows - r) / (double)rows;
		}
		output = (output + mY)*wei;
		MatrixXd temp;
		temp.setOnes(rows, cols);
		W = wei*temp;
		
	}


	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::low_rank_ssc_C_TG(MatrixXd Y, double c1, double nsig, MatrixXd &output, MatrixXd &W) {
		int rows = Y.rows();
		int cols = Y.cols();
		MatrixXd mY;
		mY.setZero(rows, cols); // mean value (normal component) of those similar patches
		for (int i = 0; i < rows; i++)
		{
			double avg = 0;
			for (int j = 0; j < cols; j++)
				avg += Y(i, j);
			avg = avg / cols;

			for (int j = 0; j < cols; j++)
				mY(i, j) = avg;
		}
		Y = Y - mY;
		JacobiSVD<MatrixXd> svdX(Y, ComputeThinU | ComputeThinV);
		MatrixXd U0 = svdX.matrixU();
		MatrixXd Sigma = svdX.singularValues();
		MatrixXd V0 = svdX.matrixV();
		MatrixXd Sigma0 = Sigma.asDiagonal();
		//MatrixXd Sigma0 = Sigma.array();

		int PatNum = cols;
		double eps = 0.00001;
		MatrixXd Temp = Sigma0.array() * Sigma0.array() - PatNum*nsig*nsig ;

		for (int i = 0; i < Temp.rows(); i++)
		{
			for (int j = 0; j < Temp.cols(); j++)
			{
				if (Temp(i, j) < 0)
					Temp(i, j) = 0;
			}

		}
		Temp = Temp.array().sqrt();
		MatrixXd thr = (c1*sqrt(PatNum)*(nsig*nsig)) / (Temp.array() + eps);
		Sigma0 = Sigma0 - thr;

		int r = 0;
		for (int i = 0; i < Sigma0.rows(); i++)
		{
			for (int j = 0; j < Sigma0.cols(); j++)
			{
				if (Sigma0(i, j) <= 0)
					Sigma0(i, j) = 0;
				else
					r = r + 1;
			}
		}
		int choice = 1; 
		if (choice == 1)
		{
			if (r == 0)
				output.setZero(rows, cols);
			else
				output = U0.middleCols(0, r - 1)*Sigma0.middleCols(0, r - 1).middleRows(0, r - 1) *V0.middleCols(0, r - 1).transpose();
		} 
		else
		{
			int number_of_out_iter = 2; //5
			int R = 10;// 20
			double rho = 1; //1
			MatrixXd X = Y.array();

			MatrixXd A, B, X_rec;
			for (int iter = 1; iter <= number_of_out_iter; iter++)
			{
				//JacobiSVD<MatrixXd> svdX(X, ComputeThinU | ComputeThinV);
				JacobiSVD<MatrixXd> svdX(X, ComputeFullU | ComputeFullV);
				MatrixXd U0 = svdX.matrixU();
				MatrixXd V0 = svdX.matrixV();
				A = U0.middleCols(0, R - 1).transpose();
				B = V0.middleCols(0, R - 1).transpose();
				X_rec = admmAXB_TG(A, B, X, rho, thr, R);
				X = X_rec.array();
			}
			output = X;
		}

		// compute weights
		double wei;
		if (r == rows)
		{
			wei = 1 / (double)rows;
		}
		else
		{
			wei = (rows - r) / (double)rows;
		}
		output = (output + mY)*wei;
		MatrixXd temp;
		temp.setOnes(rows, cols);
		W = wei*temp;

	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	MatrixXd ExKernelT<ExItems>::admmAXB_TG(MatrixXd A, MatrixXd B, MatrixXd X, double rho, MatrixXd thr, int R) {
		int MAX_ITER = 2; //5
		MatrixXd U = X.array();
		int m = X.rows();
		int n = X.cols();
		int min_m_n = m < n ? m : n;
		MatrixXd V ;
		V.setZero(m, n);
		MatrixXd Y = X.array();
		//
		MatrixXd AB = A.transpose()*B;
		double eeppss = 0.0001;
		double gamma = 0.01; //100
		double mu = 1; //1
		double gamma2 = 0.01; //0.01
		MatrixXd sigE;
		sigE.setZero(min_m_n, 1);
		for (int k = 1; k <= MAX_ITER; k++)
		{
			// U-update
			MatrixXd tem = X + V / rho;
			DC(tem, mu, thr, sigE, gamma2, U, sigE);
			// X-update
			MatrixXd lastX = X;
			for (int ii = 1; ii <= 2; ii++) //ii<=5
			{
				MatrixXd M = A*X*B.transpose();
				MatrixXd FM = (1 + gamma) / ((gamma + M.array()) + eeppss);
				
				MatrixXd temp = (A.transpose())*(FM.transpose())*B;; // AB*M.transpose();
				for (int i = 0; i < temp.rows(); i++)
				{
					for (int j = 0; j < temp.cols(); j++)
					{
						if (abs(temp(i, j)) <= 0.0001)
						{
							temp(i, j) = 0;
						}
					}
				}
				X = (mu*Y - V + temp + rho*U) / (rho + mu);
			}
			//V-update
			V = V + rho*(X - U);
			//X = U;
		}
		return X;
	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::DC(MatrixXd D, double rho, MatrixXd thr, MatrixXd T0, double a, MatrixXd &X, MatrixXd &T) {
		JacobiSVD<MatrixXd> svdX(D, ComputeThinU | ComputeThinV);
		MatrixXd Sigma = svdX.singularValues();
		MatrixXd U0 = svdX.matrixU();
		MatrixXd V0 = svdX.matrixV();
		MatrixXd thr2;
		thr2.setZero(thr.rows(), 1);
		for (int i = 0; i < thr.rows(); i++)
		{
			thr2(i, 0) = thr(i, i);
		}
		for (int k = 1; k <= 2; k++)
		{
			double lambda = rho; //0.5*rho
			MatrixXd temp = a + T.array();
			MatrixXd grad = ((1 + a)*a) / (temp.array()*temp.array() + 0.000001);
			MatrixXd T1 = Sigma.array() - lambda*(grad.array()*thr2.array());
			
			//r = 0;
			for (int i = 0; i < T1.rows(); i++)
			{
				for (int j = 0; j < T1.cols(); j++)
				{
					if (T1(i, j) < 0)
						T1(i, j) = 0;
					//else
						//r = r + 1;
				}

			}
			X = U0*T1.asDiagonal()*V0.transpose();
			double error = ((T1 - T0).array()*(T1 - T0).array()).sum();
			if (error < 0.000001)
			{
				break;
			}
			T0 = T1.array();
		}
		T = T0.array();
	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::compute_normals_NLLR() {
		//update vertex normals
		VertexIterator vi = vertex_begin();
		for (; vi != vertex_end(); ++vi)
		{
			HalfedgeHandle hh = vi->halfedge_handle_;
			VertexHandle vh = vertex_handle(vi->halfedge_handle_);
			vertex_ref(vh).normal_ = ((vertex_ref(vh).sumNormal_).normalize() / (vertex_ref(vh).weights_)).normalize();
		}

		//update face normals based on the estimated vertex normals
		FacetIterator fi = facet_begin();
		for (; fi != facet_end(); ++fi) {
			if (fi->status_.is_deleted()) continue;
			FacetHandle _fh = facet_handle(fi->halfedge_handle_);
			const HalfedgeHandle&   hh = halfedge_handle(_fh);
			const HalfedgeHandle& p_hh = prev_halfedge_handle(hh);
			const HalfedgeHandle& n_hh = next_halfedge_handle(hh);

			const Normal& cd0 = vertex_ref(vertex_handle(hh)).normal_;
			const Normal& cd1 = vertex_ref(vertex_handle(p_hh)).normal_;
			const Normal& cd2 = vertex_ref(vertex_handle(n_hh)).normal_;
			fi->normal_ = (cd0 + cd1 + cd2) / 3;
			if (fi->normal_.length() > 1E-7)
			{
				fi->normal_ = fi->normal_.normalize();
			}
			
		}
		
	}
	///////////////////////////////////////////////////////////////
	//enumerate the vertexs in counter clock fashion
	template<class ExItems>
	void ExKernelT<ExItems>::getVextexHandlesForFace(FacetHandle &_fh, std::vector<VertexHandle> &_vhs){
		std::vector<HalfedgeHandle> halfEdgeVec;	halfEdgeVec.resize(3); _vhs.resize(3);
		halfEdgeVec[0] = halfedge_handle(_fh);
		halfEdgeVec[1] = next_halfedge_handle(halfEdgeVec[0]);
		halfEdgeVec[2] = next_halfedge_handle(halfEdgeVec[1]);

		for (int i = 0; i < 3; ++i){
			VertexHandle currentVexHandle = vertex_handle(halfEdgeVec[i]);
			_vhs[i] = currentVexHandle;
		}
	}
	///////////////////////////////////////////////////////////////
	//enumerate the vertexs in counter clock fashion
	template<class ExItems>
	void ExKernelT<ExItems>::getPointsForFace(FacetHandle &_fh, std::vector<Coord> &_vhs){
		std::vector<HalfedgeHandle> halfEdgeVec;	halfEdgeVec.resize(3); _vhs.resize(3);
		halfEdgeVec[0] = halfedge_handle(_fh);
		halfEdgeVec[1] = next_halfedge_handle(halfEdgeVec[0]);
		halfEdgeVec[2] = next_halfedge_handle(halfEdgeVec[1]);

		for (int i = 0; i < 3; ++i){
			VertexHandle currentVexHandle = vertex_handle(halfEdgeVec[i]);
			_vhs[i] = coord(currentVexHandle);
		}
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getFacesForPoint(VertexHandle vh, std::vector<FacetHandle> &faceHVector){
		HalfedgeHandle& hh = halfedge_handle(vh);
		HalfedgeHandle css(hh);
		do{
			FacetHandle currentFaceHandle = facet_handle(css);
			if (currentFaceHandle.idx() == -1)
			{
				cout << "Current face handle is -1" << endl;
			}
			else{
				int ii;
				for (ii = 0; ii < faceHVector.size(); ++ii)
				{
					FacetHandle tmpFH = faceHVector[ii];
					if (currentFaceHandle.idx() == tmpFH.idx())
					{
						break;
					}
				}
				if (ii >= faceHVector.size())
				{
					faceHVector.push_back(currentFaceHandle);
				}
			}

			css = cw_rotated(css);
		} while (css != hh);
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getFaceAera(std::vector<double> &areaVector){
		//update_area();

		int faceNumber = facet_size();
		areaVector.resize(faceNumber);

		FacetIterator fit(facet_begin());
		for (; fit < facet_end(); ++fit){
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle    fh = facet_handle(hh);
			areaVector[fh.idx()] = get_area(fh);
		}
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getFaceNormal(std::vector<Normal> &normalVector){
		//update_normals();

		int faceNumber = facet_size();
		normalVector.resize(faceNumber);

		FacetIterator fit(facet_begin());
		for (; fit < facet_end(); ++fit)
		{
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle fh = facet_handle(hh);
			normalVector[fh.idx()] = normal(fh);
		}
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::UpdateVertexPosition(std::vector<Normal> &original_facet_normals, std::vector<Normal> &filtered_normals){
		update_area();
		int facet_num = facet_size();
		if (facet_num <= 0)
		{
			return;
		}
		std::vector<std::vector<Coord> > rotated_gradients;    //gradient vectors of each triangle face
		//calculate the new gradient vectors first
		calculateNewGradientVector(original_facet_normals, filtered_normals, rotated_gradients);

		//then get the divergence vector of the gradient vector
		std::vector<std::vector<double> > divergenceVector;                       //divergence vector
		calculateDivergence(rotated_gradients, divergenceVector);

		time_t start = 0, end = 0;
		//thirdly get the coefficient matrix of the Epoisson term
		start = clock();
		int vertices_number = vertex_size();
		Eigen::SparseMatrix<double> p_Coefficient_Matrix(vertices_number, vertices_number);		//the coefficient  matrix of Epoisson term
		getPoissonCoefficientMatrix(p_Coefficient_Matrix);
		end = clock();
		cout << "[Vertex updating: E_possion term] running time: " << double(end - start) / CLOCKS_PER_SEC << " second (s)" << endl;
		
		Eigen::SparseMatrix<double> total_Coefficient_Matrix = p_Coefficient_Matrix; // only with E_possion term
		Eigen::MatrixXd coordinatesVector;
		//lastly,solve the linear equation and get the result
		start = clock();
		//solveLinearEquation(total_Coefficient_Matrix, divergenceVector, coordinatesVector);
		solveLinearEquation_CG(total_Coefficient_Matrix, divergenceVector, coordinatesVector);
		end = clock();
		cout << "[Vertex updating: solve the linear equation] running time: " << double(end - start) / CLOCKS_PER_SEC << " second (s)" << endl;

		//change data type
		int verticesNum = vertex_size();
		std::vector<Coord> new_points(verticesNum);
		for (int i = 0; i < verticesNum; ++i){
			Coord iThPoint = new_points[i];
			for (int j = 0; j < 3; ++j){
				iThPoint[j] = coordinatesVector(i, j);
			}
			new_points[i] = iThPoint;
		}

		VertexIterator vit(vertex_begin());
		for (; vit != vertex_end(); ++vit){
			HalfedgeHandle currentHH = vit->halfedge_handle_;
			if (!currentHH.is_valid())
			{
				continue;
			}

			VertexHandle currentVH = vertex_handle(vit->halfedge_handle_);
			Vertex currentVex = vertex_ref(currentVH);
			currentVex.coord_ = new_points[currentVH.idx()];
			vertices_[currentVH.idx()] = currentVex;
		}

		cout << "vertex updating is ending..." << endl;
	}
	///////////////////////////////////////////////////////////////
	//calculate the new gradient vectors with the filtered normals
	template<class ExItems>
	void ExKernelT<ExItems>::calculateNewGradientVector(std::vector<Normal> &original_normals, std::vector<Normal> &filtered_normals,
		std::vector<std::vector<Coord> > &rotated_gradients){
		int facet_num = facet_size();
		rotated_gradients.resize(facet_num);

		//std::vector<Normal> original_normals;
		//getFaceNormal(original_normals);
		std::vector<double> triangle_area;
		getFaceAera(triangle_area);

		Eigen::Matrix3d rotation_matrix;
		Eigen::Vector3d vertex_vector;
		FacetIterator fit(facet_begin());
		for (; fit<facet_end(); fit++){
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle fh = facet_handle(hh);

			Normal original_normal = original_normals[fh.idx()];
			Normal filtered_normal = filtered_normals[fh.idx()];

			//get the rotation matrix for each face with original normal and filtered normal
			getRotationMatrixForEachFace(original_normal, filtered_normal, rotation_matrix);
			std::vector<Coord> points;      //store three new coordinates of the triangle
			points.resize(3);
			//obtain the three vertex coordinates
			getPointsForFace(fh, points);

			//update the three vertex coordinates of each triangle based on the original normal and the filtered normal
			for (int i = 0; i < 3; i++){
				Coord face_vertex = points[i];
				//turn the Coord type parameter to Eigen::Vector3f type
				for (int j = 0; j < 3; j++){
					vertex_vector[j] = face_vertex[j];
				}
				//convert the Eigen::Vector3f type to Coord type
				vertex_vector = rotation_matrix * vertex_vector;
				for (int j = 0; j < 3; j++){
					face_vertex[j] = vertex_vector[j];
				}
				points[i] = face_vertex;
			}

			int face_index = fh.idx();
			double area = triangle_area[face_index];
			std::vector<Coord> face_gradient_vector;       //used to store three gradient vectors for three weights
			std::vector<Coord> triangle_basis_vector;
			face_gradient_vector.resize(3);
			triangle_basis_vector.resize(3);

			for (int i = 0; i < 3; ++i){
				triangle_basis_vector[i] = ((points[(i + 1) % 3] - points[(i + 2) % 3]) % filtered_normal) / (2 * area);
			}

			//iterate the three weight of each coordinate
			for (int i = 0; i < 3; i++){
				//face_gradient_vector[i] = Vec3f(0,0,0);
				face_gradient_vector[i] = Coord(0, 0, 0);
				Coord current_vertex;
				for (int j = 0; j < 3; j++){
					current_vertex = points[j];
					face_gradient_vector[i] += current_vertex[i] * triangle_basis_vector[j];
				}
			}
			rotated_gradients[face_index] = face_gradient_vector;
		}
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getRotationMatrixForEachFace(Normal &original_normal,
		Normal &filtered_normal, Eigen::Matrix3d &rotation_matrix){
		Normal unit_original_normal = original_normal.normalize();
		Normal unit_filtered_normal = filtered_normal.normalize();

		Coord unit_rotation_axis = unit_original_normal % unit_filtered_normal;
		double sine_value = unit_rotation_axis.length();
		double cosine_value = unit_original_normal * unit_filtered_normal;

		rotation_matrix.resize(3, 3);
		rotation_matrix(0, 0) = 0;
		rotation_matrix(0, 1) = -unit_rotation_axis[2];
		rotation_matrix(0, 2) = unit_rotation_axis[1];
		rotation_matrix(1, 0) = unit_rotation_axis[2];
		rotation_matrix(1, 1) = 0;
		rotation_matrix(1, 2) = -unit_rotation_axis[0];
		rotation_matrix(2, 0) = -unit_rotation_axis[1];
		rotation_matrix(2, 1) = unit_rotation_axis[0];
		rotation_matrix(2, 2) = 0;

		rotation_matrix = rotation_matrix + (rotation_matrix * rotation_matrix) / (1 + cosine_value);
		rotation_matrix(0, 0) += 1;
		rotation_matrix(1, 1) += 1;
		rotation_matrix(2, 2) += 1;
	}

	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::calculateDivergence(std::vector<std::vector<Coord> >&rotated_gradients,
		std::vector<std::vector<double> > &divergenceVector){
		int n_vertices = vertex_size();
		divergenceVector.resize(n_vertices);
		for (size_t i = 0; i < divergenceVector.size(); ++i){
			std::vector<double> zeroVector(3, 0);
			divergenceVector[i] = zeroVector;
		}

		std::vector<Normal> normal_vector;
		getFaceNormal(normal_vector);

		std::vector<double> face_areas;
		getFaceAera(face_areas);

		std::vector<Coord> current_gradient_vec;
		current_gradient_vec.resize(3);

		double areaSum = 0.0;
		for (int i = 0; i < (int)face_areas.size(); ++i){
			areaSum += face_areas[i];
		}
		double averageArea = areaSum / (int)facet_size();

		FacetIterator fit(facet_begin());
		for (; fit<facet_end(); fit++){
			//È¡³öÃ¿¸öÃæµÄ¾ä±ú
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle fh = facet_handle(hh);

			std::vector<Coord> points;      //store three new coordinates of the triangle
			std::vector<int> vertexIndexVector;
			points.resize(3);   vertexIndexVector.resize(3);
			int index = 0;

			std::vector<VertexHandle> vertexHandlesVector;	vertexHandlesVector.resize(3);
			getVextexHandlesForFace(fh, vertexHandlesVector);
			for (size_t i = 0; i < vertexHandlesVector.size(); ++i)
			{
				VertexHandle currentVH = vertexHandlesVector[i];
				points[i] = coord(currentVH);
				vertexIndexVector[i] = currentVH.idx();
			}

			int face_index = fh.idx();
			double area = face_areas[face_index];
			std::vector<Coord> triangle_basis_vector;
			triangle_basis_vector.resize(3);

			Normal currentNormal = normal_vector[face_index];

			for (int i = 0; i < 3; ++i){
				triangle_basis_vector[i] = ((points[(i + 1) % 3] - points[(i + 2) % 3]) % currentNormal) / (2 * area);
			}

			std::vector<Coord> currentNewGradients = rotated_gradients[face_index];

			for (int i = 0; i < 3; i++){
				Coord currentBasis;
				int currentVertexIndex;
				Coord currentGradient = currentNewGradients[i];

				for (int j = 0; j < 3; j++){
					currentBasis = triangle_basis_vector[j];
					currentVertexIndex = vertexIndexVector[j];

					std::vector<double> currentDivergence = divergenceVector[currentVertexIndex];
					double dotVector = (currentGradient * currentBasis);
					currentDivergence[i] += dotVector * area / averageArea;
					divergenceVector[currentVertexIndex] = currentDivergence;
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getPoissonCoefficientMatrix(Eigen::SparseMatrix<double> &p_Coefficient_Matrix){
		//resize the coefficient matrix size and then set it zero
		double vertices_number = vertex_size();
		p_Coefficient_Matrix.resize(vertices_number, vertices_number);

		std::vector<double> face_area(facet_size());
		getFaceAera(face_area);

		//calculate the average area of the mesh
		double areaSum = 0.0;
		for (int i = 0; i < (int)face_area.size(); ++i){
			areaSum += face_area[i];
		}
		double averageArea = areaSum / facet_size();

		time_t start = 0, end = 0; //lzhu
		time(&start);//lzhu
		FacetIterator fit(facet_begin());
		for (; fit<facet_end(); fit++){
			
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle    fh = facet_handle(hh);

			int index = 0;
			std::vector<Coord> pointVector;    std::vector<int> vertexIndexVector;
			pointVector.resize(3);  vertexIndexVector.resize(3);

			std::vector<VertexHandle> vertexHandleVector;
			getVextexHandlesForFace(fh, vertexHandleVector);
			for (size_t i = 0; i < vertexHandleVector.size(); ++i)
			{
				VertexHandle currentVH = vertexHandleVector[i];
				pointVector[i] = coord(currentVH);
				vertexIndexVector[i] = currentVH.idx();
			}

			std::vector<Coord> edgeVector;  edgeVector.resize(3);
			for (int i = 0; i < 3; ++i){
				edgeVector[i] = pointVector[(i + 2) % 3] - pointVector[(i + 1) % 3];
			}

			double currentArea = face_area[fh.idx()];
			for (int i = 0; i < 3; ++i){
				int firstPointIndex = vertexIndexVector[i];
				Coord firstEdge = edgeVector[i];
				for (int j = 0; j < 3; ++j){
					int secondPointIndex = vertexIndexVector[j];
					Coord secondEdge = edgeVector[j];
					double dotVector = firstEdge * secondEdge;
					//               p_Coefficient_Matrix.coeffRef(firstPointIndex,secondPointIndex) += dotVector / (4 * currentArea);
					p_Coefficient_Matrix.coeffRef(firstPointIndex, secondPointIndex) += dotVector / (4 * currentArea * averageArea);
				}
			}
		}
		time(&end);
		cout << " [Constructing poisson coefficient matrix]running time: " << (end - start) << " second (s)" << endl; //lzhu
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::calculateFaceRotationMatrix(Coord rotationAxisVec,
		Normal firstFaceNormal, Normal secondFaceNormal, Eigen::Matrix3d &rotationMatrix){
		double cosine_theta = (firstFaceNormal * secondFaceNormal) / (firstFaceNormal.norm() * secondFaceNormal.norm());
		double sine_theta = sqrt(1 - cosine_theta * cosine_theta);

		//TriMesh::Point rotationAxisVector = rotationAxisVec.normalize();
		Coord rotationAxisVector = rotationAxisVec.normalize();
		double x = rotationAxisVector[0], y = rotationAxisVector[1], z = rotationAxisVector[2];

		rotationMatrix(0, 0) = x * x + (1 - x * x) * cosine_theta;
		rotationMatrix(0, 1) = x * y * (1 - cosine_theta) - z * sine_theta;
		rotationMatrix(0, 2) = x * z * (1 - cosine_theta) + y * sine_theta;
		rotationMatrix(1, 0) = y * x * (1 - cosine_theta) + z * sine_theta;
		rotationMatrix(1, 1) = y * y + (1 - y * y) * cosine_theta;
		rotationMatrix(1, 2) = y * z * (1 - cosine_theta) - x * sine_theta;
		rotationMatrix(2, 0) = z * x * (1 - cosine_theta) - y * sine_theta;
		rotationMatrix(2, 1) = z * y * (1 - cosine_theta) + x * sine_theta;
		rotationMatrix(2, 2) = z * z + (1 - z * z) * cosine_theta;
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::updateCoefficientMatrixWithQuadraticTerm(int faceHandleIndex,
		double edgeLength, double weightedValue, Eigen::SparseMatrix<double> &s_Coefficient_Matrix){
		//TriMesh::FaceHandle currentFaceHandle = mesh.face_handle(faceHandleIndex);
		FacetHandle currentFaceHandle = FacetHandle(faceHandleIndex);
		std::vector<Coord> pointVector;    std::vector<int> vertexIndexVector;
		pointVector.resize(3);		vertexIndexVector.resize(3);
		int index = 0;

		std::vector<VertexHandle> vhVector;
		getVextexHandlesForFace(currentFaceHandle, vhVector);
		for (size_t i = 0; i < vhVector.size(); ++i)
		{
			VertexHandle tmpVH = vhVector[i];
			pointVector[index] = coord(tmpVH);
			vertexIndexVector[index] = tmpVH.idx();
			index++;
		}
		//enumerate vertexs of the face in counter clockwise fashion
		std::vector<Coord> edgeVector; edgeVector.resize(3);
		edgeVector[0] = pointVector[2] - pointVector[1];
		edgeVector[1] = pointVector[0] - pointVector[2];
		edgeVector[2] = pointVector[1] - pointVector[0];

		double triangleArea = 0.5 * (edgeVector[0] % edgeVector[1]).length();

		for (int i = 0; i < 3; i++){
			int firstVertexIndex = vertexIndexVector[i];
			for (int j = 0; j < 3; j++){
				int secondVertexIndex = vertexIndexVector[j];
				s_Coefficient_Matrix.coeffRef(firstVertexIndex, secondVertexIndex) += (edgeVector[i] * edgeVector[j]) * edgeLength * weightedValue / (4 * triangleArea * triangleArea);
				//s_Coefficient_Matrix(firstVertexIndex,secondVertexIndex) += (edgeVector[i] * edgeVector[j]) * edgeLength / (4 * triangleArea);
			}
		}
	}

	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::updateCoefficientMatrixWithCrossTerm(int firstFaceHandleIndex,
		int secondFaceHandleIndex, int conjunctEdgeIndex, double weightedValue, Eigen::SparseMatrix<double> &s_Coefficient_Matrix){
		FacetHandle firstFaceHandle = FacetHandle(firstFaceHandleIndex);
		FacetHandle secondFaceHandle = FacetHandle(secondFaceHandleIndex);
		EdgeHandle conjunctEdgeHandle = EdgeHandle(conjunctEdgeIndex);

		VertexHandle firstConjunctPointHandle = vertex_handle(halfedge_handle(conjunctEdgeHandle, 0));
		VertexHandle secondConjunctPointHandle = vertex_handle(halfedge_handle(conjunctEdgeHandle, 1));

		Coord firstConjunctPoint = coord(firstConjunctPointHandle);
		Coord secondConjunctPoint = coord(secondConjunctPointHandle);
		Coord rotationAxisVec = firstConjunctPoint - secondConjunctPoint;

		Normal firstFaceNormal = normal(firstFaceHandle);
		Normal secondFaceNormal = normal(secondFaceHandle);

		Eigen::Matrix3d rotationMatrix;
		calculateFaceRotationMatrix(rotationAxisVec, firstFaceNormal, secondFaceNormal, rotationMatrix);

		//used to store the three coordinates of two triangle faces that are linked with the conjunct edge
		std::vector<Coord> firstFacePointVec, secondFacePointVec;
		firstFacePointVec.resize(3); secondFacePointVec.resize(3);
		std::vector<int> firstFacePointIndexVec, secondFacePointIndexVec;
		firstFacePointIndexVec.resize(3); secondFacePointIndexVec.resize(3);

		int index = 0;
		std::vector<VertexHandle> firstFaceVHs;
		getVextexHandlesForFace(firstFaceHandle, firstFaceVHs);
		for (size_t i = 0; i < firstFaceVHs.size(); ++i)
		{
			VertexHandle currentHandle = firstFaceVHs[i];
			firstFacePointIndexVec[index] = currentHandle.idx();
			firstFacePointVec[index] = coord(currentHandle);
			index++;
		}
		index = 0;
		std::vector<VertexHandle> secondFaceVHs;
		getVextexHandlesForFace(secondFaceHandle, secondFaceVHs);
		for (size_t j = 0; j < secondFaceVHs.size(); ++j)
		{
			VertexHandle currentHandle = secondFaceVHs[j];
			secondFacePointIndexVec[index] = currentHandle.idx();
			secondFacePointVec[index] = coord(currentHandle);
			index++;
		}
		//TriMesh::Point type vector
		std::vector<Coord> firstGradientWeightVec, secondGradientWeightVec;
		firstGradientWeightVec.resize(3); secondGradientWeightVec.resize(3);

		for (int i = 0; i < 3; i++){
			firstGradientWeightVec[i] = firstFacePointVec[(i + 2) % 3] - firstFacePointVec[(i + 1) % 3];
			secondGradientWeightVec[i] = secondFacePointVec[(i + 2) % 3] - secondFacePointVec[(i + 1) % 3];
		}

		double firstFaceArea = 0.5 * (firstGradientWeightVec[0] % firstGradientWeightVec[1]).length();
		double secondFaceArea = 0.5 * (secondGradientWeightVec[0] % secondGradientWeightVec[1]).length();
		double conjunctEdgeLength = rotationAxisVec.length();

		for (int i = 0; i < 3; ++i){
			firstGradientWeightVec[i] = (firstFaceNormal % firstGradientWeightVec[i]) * (double)(0.5 / firstFaceArea);
			secondGradientWeightVec[i] = (secondFaceNormal % secondGradientWeightVec[i]) * (double)(0.5 / secondFaceArea);
		}

		//change the form of gradient weight vectors
		std::vector<Eigen::Vector3d> firstEigenGradientWeightVec, secondEigenGradientWeightVec;
		firstEigenGradientWeightVec.resize(3); secondEigenGradientWeightVec.resize(3);
		for (int i = 0; i < 3; ++i){
			Eigen::Vector3d firstEigenGradientWeight = firstEigenGradientWeightVec[i], secondEigenGradientWeight = secondEigenGradientWeightVec[i];
			Coord firstPointGradientWeight = firstGradientWeightVec[i], secondPointGradientWeight = secondGradientWeightVec[i];
			for (int j = 0; j < 3; ++j){
				firstEigenGradientWeight(j) = firstPointGradientWeight[j];
				secondEigenGradientWeight(j) = secondPointGradientWeight[j];
			}
		}

		for (int i = 0; i < 3; ++i){
			int rowIndex = firstFacePointIndexVec[i];
			Eigen::Vector3d rowWeightVec = firstEigenGradientWeightVec[i];
			for (int j = 0; j < 3; j++){
				int columnIndex = secondFacePointIndexVec[j];
				Eigen::Vector3d columnWeightVec = secondEigenGradientWeightVec[j];
				s_Coefficient_Matrix.coeffRef(rowIndex, columnIndex) += (-2 * conjunctEdgeLength) * weightedValue * rowWeightVec.transpose() * rotationMatrix * columnWeightVec;
				//s_Coefficient_Matrix(rowIndex,columnIndex) += (-2 * conjunctEdgeLength) * rowWeightVec.transpose() * rotationMatrix * columnWeightVec;
			}
		}
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getSmoothCoefficientMatrix(Eigen::SparseMatrix<double> &s_Coefficient_Matrix){
		int vertices_number = vertex_size();
		s_Coefficient_Matrix.resize(vertices_number, vertices_number);
		//s_Coefficient_Matrix.setZero();

		double averageEdgeLength = get_average_edge_length(), lanmada = 0.5;

		EdgeIterator eit = edge_begin();
		for (; eit != edge_end(); ++eit) {
			Halfedge firstEdge = eit->halfedges_[0];
			VertexHandle tmpVH = firstEdge.vertex_handle_;		//the to point of the first edge
			HalfedgeHandle tmpHalfEH = halfedge_handle(tmpVH);		//the half edge of the firstEdge represent
			EdgeHandle currentEdgeHandle = edge_handle(tmpHalfEH);

			//cout << "smooth coefficient matrix --- edge index : " << currentEdgeHandle.idx() << endl;
			//if (is_boundary_edge(currentEdgeHandle))
			if (!is_boundary_edge(currentEdgeHandle)) //lzhu
			{
				HalfedgeHandle firstHalfEdgeHandle = halfedge_handle(currentEdgeHandle, 0);
				HalfedgeHandle secondHalfEdgeHandle = halfedge_handle(currentEdgeHandle, 1);

				FacetHandle firstFaceHandle = facet_handle(firstHalfEdgeHandle);
				FacetHandle secondFaceHandle = facet_handle(secondHalfEdgeHandle);

				if (firstFaceHandle.idx() == -1 || secondFaceHandle.idx() == -1)
				{
					continue;
				}

				VertexHandle firstPointHandle = vertex_handle(firstHalfEdgeHandle);
				VertexHandle secondPointHandle = vertex_handle(secondHalfEdgeHandle);

				Coord firstPointCoord = coord(firstPointHandle);
				Coord secondPointCoord = coord(secondPointHandle);

				updateCoefficientMatrixWithQuadraticTerm(firstFaceHandle.idx(),
					(firstPointCoord - secondPointCoord).length(), (lanmada / averageEdgeLength), s_Coefficient_Matrix);
				updateCoefficientMatrixWithQuadraticTerm(secondFaceHandle.idx(),
					(firstPointCoord - secondPointCoord).length(), (lanmada / averageEdgeLength), s_Coefficient_Matrix);
				updateCoefficientMatrixWithCrossTerm(firstFaceHandle.idx(), secondFaceHandle.idx(),
					currentEdgeHandle.idx(), (lanmada / averageEdgeLength), s_Coefficient_Matrix);
			}
		}
	}
	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::solveLinearEquation(Eigen::SparseMatrix<double> &t_Coefficient_Matrix,
		std::vector<std::vector<double> > &divergenceVector, Eigen::MatrixXd &coordinatesVector){
		int rowNum = t_Coefficient_Matrix.rows(), colNum = t_Coefficient_Matrix.cols();
		if (rowNum != colNum)
			return;
		//turn the std::vector type data into Eigen::MatrixXd type
		int vectorSize = divergenceVector.size();
		if (rowNum != vectorSize)
			return;

		Eigen::MatrixXd divergenceMatrix(vectorSize, 3);
		for (int i = 0; i < vectorSize; ++i){
			std::vector<double> iThDivergenceVec = divergenceVector[i];
			for (int j = 0; j < 3; ++j){
				divergenceMatrix(i, j) = iThDivergenceVec[j];
			}
		}

		//Eigen::MatrixXd coordinatesVector(vectorSize,3);
		//coordinatesVector = t_Coefficient_Matrix.llt().solve(divergenceMatrix);

		coordinatesVector.resize(vectorSize, 3);
		//cout << "rows=" << coordinatesVector.rows() << ", cols=" << coordinatesVector.cols() << endl;

		//coordinatesVector = t_Coefficient_Matrix.colPivHouseholderQr().solve(divergenceVector);

		// solve linear system
		Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
		solver.analyzePattern(t_Coefficient_Matrix);
		solver.factorize(t_Coefficient_Matrix);
		coordinatesVector = solver.solve(divergenceMatrix);
	}

	///////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::solveLinearEquation_CG(Eigen::SparseMatrix<double> &t_Coefficient_Matrix,
		std::vector<std::vector<double> > &divergenceVector, Eigen::MatrixXd &coordinatesVector){
		int rowNum = t_Coefficient_Matrix.rows(), colNum = t_Coefficient_Matrix.cols();
		if (rowNum != colNum)
			return;
		//turn the std::vector type data into Eigen::MatrixXd type
		int vectorSize = divergenceVector.size();
		if (rowNum != vectorSize)
			return;

		Eigen::MatrixXd divergenceMatrix(vectorSize, 3);
		for (int i = 0; i < vectorSize; ++i){
			std::vector<double> iThDivergenceVec = divergenceVector[i];
			for (int j = 0; j < 3; ++j){
				divergenceMatrix(i, j) = iThDivergenceVec[j];
			}
		}

		//Eigen::MatrixXd coordinatesVector(vectorSize,3);
		//coordinatesVector = t_Coefficient_Matrix.llt().solve(divergenceMatrix);
		ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
		coordinatesVector = solver.compute(t_Coefficient_Matrix).solve(divergenceMatrix);
	}


	////////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::adjustVertexCoord(int iterations){///Sun Xiangfang TVCG 2007 fast and effective feature-preserving mesh denoising			
			int vertex_num = vertex_size();
			while(--iterations)
			{	
				std::vector<Coord> updateVertexPosition;
				updateVertexPosition.resize(vertex_num);
				for (int i = 0; i < vertex_num; i++){
					//cout << i << endl;
					VertexHandle vh(i);
					Coord&       vc = coord(vh);
					HalfedgeHandle& hh = halfedge_handle(vh);
					HalfedgeHandle  css(hh);
					do{
						HalfedgeHandle opp_hh = opposite_halfedge_handle(css);
						Coord&         opp_vc = coord(vertex_handle(opp_hh));
						FacetHandle    fl = facet_handle(css);
						FacetHandle    fr = facet_handle(opp_hh);

						if (fl.is_valid()){
							updateVertexPosition[i] += facet_ref(fl).normal_*(facet_ref(fl).normal_ *(opp_vc - vc));

						}
						if (fr.is_valid()){
							updateVertexPosition[i] += facet_ref(fr).normal_*(facet_ref(fr).normal_ *(opp_vc - vc));

						}

						css = cw_rotated(css);

					} while (css != hh);
				}
				for (int i = 0; i < vertex_num; i++)
				{ 
					vertex_ref(VertexHandle(i)).coord_ += updateVertexPosition[i] * 1.0 / 18.0; 
				}
			};
		}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::calcUpdatedVertexNormal(std::vector<Normal>& updatedVertexNormals)
	{
			VertexIterator vi = vertex_begin();
			for (; vi != vertex_end(); vi++)
			{
				Normal updatedNormal = vi->normal_;
				updatedVertexNormals.push_back(updatedNormal);
			}
	}
	
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Scalar
		ExKernelT<ExItems>::calc_facet_area(const FacetHandle& _fh){

			HalfedgeHandle& hh0 = halfedge_handle(_fh);
			HalfedgeHandle& hh1 = next_halfedge_handle(hh0);
			HalfedgeHandle& hh2 = next_halfedge_handle(hh1);

			EdgeHandle& eh0 = edge_handle(hh0);
			EdgeHandle& eh1 = edge_handle(hh1);
			EdgeHandle& eh2 = edge_handle(hh2);

			Scalar eh0_length = calc_edge_length(eh0);
			Scalar eh1_length = calc_edge_length(eh1);
			Scalar eh2_length = calc_edge_length(eh2);

			Scalar area;
			Scalar p = (eh0_length + eh1_length + eh2_length)/2.0;
			area=sqrt( p * (p-eh0_length) * (p-eh1_length) * (p-eh2_length) );

			facet_ref(_fh).area_ = area;
			return area;
	}

	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::update_area()
	{   
		set_isArea(true);
		FacetIterator fit(facet_begin() );
		for(;fit<facet_end(); fit++)
		{
			HalfedgeHandle hh = fit->halfedge_handle_;
			FacetHandle    fh = facet_handle(hh);
			calc_facet_area(fh);
		}
	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Scalar
		ExKernelT<ExItems>::get_area(const FacetHandle& _fh)
	{
		return facet_ref(_fh).area_;
	}
	//////////////////////////////////////////////////////////////////////////////
	template<class ExItems> 
	double 
		ExKernelT<ExItems>::calc_edge_length(const EdgeHandle &_eh) {//updating all edge length


			HalfedgeHandle& h1 = halfedge_handle(_eh,0);
			HalfedgeHandle& h2 = halfedge_handle(_eh,1);

			Vertex v0 = vertex_ref(vertex_handle(h1) );
			Vertex v1 = vertex_ref(vertex_handle(h2) );

			return (v0.coord_-v1.coord_).norm();

	}

	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getNeighborVertices_withRingIndex(VertexHandle& _vh, int _verticesNum, std::vector<VertexHandle>& NeighborVertices, std::vector<int>& ringIdx)
	{
		int ringIndex = 0;
		NeighborVertices.push_back(_vh);
		ringIdx.push_back(ringIndex);
		int iteration = 0;
		int verNewNum = NeighborVertices.size();
		int verOldNum = verNewNum - 1;
		int verOldNum1 = verOldNum;

		while (1)
		{
			ringIndex += 1;
			verOldNum = NeighborVertices.size();
			for (int i = verOldNum1; i < verNewNum; i++) {
				VertexHandle& vh = NeighborVertices[i];
				HalfedgeHandle& hh = halfedge_handle(vh);
				HalfedgeHandle css(hh);
				do {
					int ii = 0;
					HalfedgeHandle& opp_hh = opposite_halfedge_handle(css);
					VertexHandle&   opp_vh = vertex_handle(opp_hh);
					for (ii = 0; ii<NeighborVertices.size(); ii++)
						if (opp_vh == NeighborVertices[ii]) break;

					if (ii >= NeighborVertices.size())
					{
						NeighborVertices.push_back(opp_vh);
						ringIdx.push_back(ringIndex);
						if (NeighborVertices.size() == _verticesNum)
						{
							return;
						}
					}

					css = cw_rotated(css);
				} while (css != hh);
			}
			verNewNum = NeighborVertices.size();
			verOldNum1 = verOldNum;
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::getNeighborVertices(VertexHandle& _vh, int _verticesNum, std::vector<VertexHandle>& NeighborVertices) {

			NeighborVertices.push_back(_vh);
			int iteration = 0;
			int verNewNum = NeighborVertices.size();
			int verOldNum = verNewNum - 1;
			int verOldNum1 = verOldNum;

			while (1)
			{
				verOldNum = NeighborVertices.size();
				for (int i = verOldNum1; i < verNewNum; i++) {
					VertexHandle& vh = NeighborVertices[i];
					HalfedgeHandle& hh = halfedge_handle(vh);
					HalfedgeHandle css(hh);
					do {
						int ii = 0;
						HalfedgeHandle& opp_hh = opposite_halfedge_handle(css);
						VertexHandle&   opp_vh = vertex_handle(opp_hh);
						for (ii = 0; ii<NeighborVertices.size(); ii++)
						if (opp_vh == NeighborVertices[ii]) break;

						if (ii >= NeighborVertices.size())
						{
							NeighborVertices.push_back(opp_vh);
							if (NeighborVertices.size() >= _verticesNum)
							{
								return;
							}
						}
						if (NeighborVertices.size() > _verticesNum) break;

						css = cw_rotated(css);
					} while (css != hh);
					if (NeighborVertices.size() > _verticesNum) break;
				}
				if (NeighborVertices.size() > _verticesNum) break;
				verNewNum = NeighborVertices.size();
				verOldNum1 = verOldNum;
			}

		}
	/////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::getNeighborRing(VertexHandle& _vh, int _ring, std::vector<VertexHandle>& NeighborRing){//µÃµ½µÚring»·ÁÚ½Óµã

			NeighborRing.push_back( _vh );
			int iteration = 0;
			int verNewNum = NeighborRing.size();
			int verOldNum = verNewNum-1;
			int verOldNum1 = verOldNum;

			do{
				verOldNum = NeighborRing.size();
				for(int i=verOldNum1; i<verNewNum; i++){
					VertexHandle& vh = NeighborRing[i];					
					HalfedgeHandle& hh = halfedge_handle(vh);
					HalfedgeHandle css(hh);
					do{
						int ii = 0;
						HalfedgeHandle& opp_hh = opposite_halfedge_handle(css);
						VertexHandle&   opp_vh = vertex_handle(opp_hh);
						for(ii=0; ii<NeighborRing.size(); ii++)
							if(opp_vh == NeighborRing[ii] ) break;

						if(ii >= NeighborRing.size() )
							NeighborRing.push_back(opp_vh);

						css = cw_rotated(css);
					}while(css != hh);
				}

				verNewNum = NeighborRing.size();
				verOldNum1 = verOldNum;
				iteration++;
			} while (iteration < _ring);
			//
	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::getNeighborFaceN1(FacetHandle& _fh, std::vector<FacetHandle>& _fhs){//sharing common edges

			HalfedgeHandle& hh = halfedge_handle(_fh);
			HalfedgeHandle& pre_hh = prev_halfedge_handle(hh);
			HalfedgeHandle& nex_hh = next_halfedge_handle(hh);

			_fhs.push_back( facet_handle( opposite_halfedge_handle(hh) ) );
			_fhs.push_back( facet_handle( opposite_halfedge_handle(pre_hh ) ) );
			_fhs.push_back( facet_handle(opposite_halfedge_handle(nex_hh ) ) );
	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void 
		ExKernelT<ExItems>::getNeighborFaceN2(FacetHandle& _fh, std::vector<FacetHandle>& _fhs){//sharing common vertices

			HalfedgeHandle& hh = halfedge_handle(_fh);
			HalfedgeHandle& pre_hh = prev_halfedge_handle(hh);
			HalfedgeHandle& nex_hh = next_halfedge_handle(hh);

			VertexHandle  vhs[3];
			vhs[0] = vertex_handle(hh);
			vhs[1] = vertex_handle(pre_hh);
			vhs[2] = vertex_handle(nex_hh);

			for(int i=0; i<3; i++){

				HalfedgeHandle& hhv = halfedge_handle( vhs[i] );
				HalfedgeHandle cursor(hhv);

				do{

					FacetHandle fh = facet_handle(cursor);
					if(fh.is_valid() && fh != _fh){

						if(_fhs.size() != 0){

							for(int j=0; j< _fhs.size(); j++){

								if(fh.idx() == _fhs[j].idx() ) break;
							}//end for

							if(j>= _fhs.size() ) _fhs.push_back(fh);

						}//end if
						else _fhs.push_back(fh);
					}//end if

					cursor = cw_rotated(cursor);

				}while(hhv != cursor);//end for do while
			}//end for

	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::output_to_file(){

		FILE *fp;
		fp=fopen("result.off","w");

		int no_vertex=vertex_size();
		int no_facet=facet_size();
		int edge = 0;

		fprintf(fp,"OFF\n");
		fprintf(fp,"%d  %d %d\n",no_vertex,no_facet, edge);

		VertexIterator vit(vertex_begin());

		for(;vit!=vertex_end();vit++){

			VertexHandle vh=vertex_handle(vit->halfedge_handle_);
			fprintf(fp," %f  %f  %f\n",coord(vh).data_[0],coord(vh).data_[1],coord(vh).data_[2]);

		}

		FacetIterator fit(facet_begin());

		for(;fit!=facet_end();fit++){

			HalfedgeHandle hh=fit->halfedge_handle_;
			HalfedgeHandle nh=next_halfedge_handle(hh);
			HalfedgeHandle nnh=next_halfedge_handle(nh);

			VertexHandle vh=vertex_handle(hh);
			VertexHandle nvh=vertex_handle(nh);
			VertexHandle nnvh=vertex_handle(nnh);

			fprintf(fp,"3 %d  %d  %d\n",vh.idx(),nvh.idx(),nnvh.idx());

		}

		fclose(fp);
	}
	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void ExKernelT<ExItems>::output_to_file(char* filename){

		FILE *fp;
		fp = fopen(filename, "w");

		int no_vertex = vertex_size();
		int no_facet = facet_size();
		int edge = 0;

		fprintf(fp, "OFF\n");
		fprintf(fp, "%d  %d %d\n", no_vertex, no_facet, edge);

		VertexIterator vit(vertex_begin());

		for (; vit != vertex_end(); vit++){

			VertexHandle vh = vertex_handle(vit->halfedge_handle_);
			fprintf(fp, " %f  %f  %f\n", coord(vh).data_[0], coord(vh).data_[1], coord(vh).data_[2]);

		}

		FacetIterator fit(facet_begin());

		for (; fit != facet_end(); fit++){

			HalfedgeHandle hh = fit->halfedge_handle_;
			HalfedgeHandle nh = next_halfedge_handle(hh);
			HalfedgeHandle nnh = next_halfedge_handle(nh);

			VertexHandle vh = vertex_handle(hh);
			VertexHandle nvh = vertex_handle(nh);
			VertexHandle nnvh = vertex_handle(nnh);

			fprintf(fp, "3 %d  %d  %d\n", vh.idx(), nvh.idx(), nnvh.idx());

		}

		fclose(fp);
	}
//////////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::normal(const FacetHandle& _fh) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );
			return facet_ref(_fh).normal_;
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::calc_normal(const FacetHandle& _fh) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );

			const HalfedgeHandle&   hh = halfedge_handle(_fh);
			const HalfedgeHandle& p_hh = prev_halfedge_handle(hh);
			const HalfedgeHandle& n_hh = next_halfedge_handle(hh);

			const Coord& cd0 = coord( vertex_handle( hh) );
			const Coord& cd1 = coord( vertex_handle(p_hh) );
			const Coord& cd2 = coord( vertex_handle(n_hh) );

			//return ((cd1-cd0)%(cd2-cd1)).normalize();
			return ((cd2-cd1)%(cd1-cd0)).normalize();//be careful
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_facet_normals(void) {

		set_isNormal(true);
		FacetIterator fi = facet_begin();

		for ( ; fi!=facet_end(); ++fi) {
			if (fi->status_.is_deleted()) continue;

			assert(fi->halfedge_handle_.is_valid());
			fi->normal_ = calc_normal( facet_handle(fi->halfedge_handle_) );
			fi->oriNormal_ = fi->normal_;
		}
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::normal(const VertexHandle& _vh) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );
			return vertex_ref(_vh).normal_;
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::calc_normal(const VertexHandle& _vh) {
			assert( _vh.is_valid());
			assert( _vh.idx() < vertex_size() );

			Normal          norm;
			HalfedgeHandle& hh = halfedge_handle(_vh);
			HalfedgeHandle  cs (hh);
			do {
				FacetHandle& fh = facet_handle(cs);
				if (fh.is_valid()) 
					norm  += normal(fh);

				cs = cw_rotated(cs);
			} while ( cs != hh );

			return norm.normalize();
	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_vertex_normals(void) {
		VertexIterator vi = vertex_begin();

		for ( ; vi!=vertex_end(); ++vi) {
			if (vi->status_.is_deleted()) continue;

			assert(vi->halfedge_handle_.is_valid());
			vi->normal_ = calc_normal( vertex_handle(vi->halfedge_handle_) );
			//std::cout<<vi->normal_<<std::endl;
		}
	}


	///////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Normal
		ExKernelT<ExItems>::calc_normal_max(const VertexHandle& _vh){

			assert(_vh.is_valid() );
			assert(_vh.idx() < vertex_size() );

			HalfedgeHandle& hh = halfedge_handle(_vh);
			Coord& vc = vertex_ref(_vh).coord_;
			Normal n(0,0,0);
			HalfedgeHandle css(hh);
			do{
				FacetHandle& fh = facet_handle(css);
				if(fh.is_valid() ){

					HalfedgeHandle& nexhh = next_halfedge_handle(css);
					HalfedgeHandle& prehh = prev_halfedge_handle(css);

					VertexHandle& nvh = vertex_handle(nexhh);
					VertexHandle& prevh = vertex_handle(prehh);

					Coord& nvc = vertex_ref(nvh).coord_;
					Coord& prevc = vertex_ref(prevh).coord_;

					Coord vec1 = vc - nvc;
					Coord vec2 = vc - prevc;
					Coord vec12cross = vec1 % vec2;//cross multiplication

					float weight = vec12cross.length() / (vec1.sqLength() * vec2.sqLength() );

					n += facet_ref(fh).normal_ * weight;//---------------------------------------×¢Òâ
					//n += calc_normal(fh);//---------------------------------×¢Òâ

				}//if

				css = cw_rotated(css);
			}while(css != hh);

			n.normalize();

			return n;


	}
	//////////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_vertex_normals_max(void){

		VertexIterator vi = vertex_begin();

		for ( ; vi!=vertex_end(); ++vi) {
			if (vi->status_.is_deleted() ) continue;

			assert(vi->halfedge_handle_.is_valid());
			//vi->normal_ = calc_normal( vertex_handle(vi->halfedge_handle_));
			vi->normal_ = calc_normal_max( vertex_handle(vi->halfedge_handle_) );
			vi->oriNormal_ = vi->normal_; //lzhu 
			//std::cout<<vi->normal_<<std::endl;
		}
	}
	///////////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_normals(void) {
		update_facet_normals();

		VertexIterator vi = vertex_begin();
		for ( ; vi!=vertex_end(); ++vi) {
			if (vi->status_.is_deleted()) continue;

			assert(vi->halfedge_handle_.is_valid());

			HalfedgeHandle& hh = vi->halfedge_handle_;
			VertexHandle&   vh = vertex_handle(hh);
			HalfedgeHandle  cs (hh);
			Normal          norm;
			do {
				const FacetHandle& fh = facet_handle(cs);
				if (fh.is_valid())  
					norm +=  facet_ref(fh).normal_;
				cs = cw_rotated( cs );
			} while ( cs != hh ); 

			vi->normal_ = norm.normalize();
		}


	}


	///////////////////////////////////////////////////////////////////////////////
	template <class ExItems> 
	void ExKernelT<ExItems>::update_edge_length(void) {//¼ÆËãÍø¸ñÖÐ±ß³¤ÐÅÏ¢
		float global_max_edge_length_ = 0;
		float averagedlength = 0.0;

		EdgeIterator eit = edge_begin();
		for ( ; eit!=edge_end(); ++eit ) {
			Vertex& v0 = vertex_ref( eit->halfedges_[0].vertex_handle_ );
			Vertex& v1 = vertex_ref( eit->halfedges_[1].vertex_handle_ );

			eit->length_ = (v0.coord_ - v1.coord_).norm();
			//std::cout<<eit->length_<<"\n";
			averagedlength += eit->length_;

			if (global_max_edge_length_ < eit->length_)
				global_max_edge_length_ = eit->length_; 
		} 

		averagedlength /= edge_size();
		set_average_edge_length(averagedlength);
	}
	/////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	typename ExKernelT<ExItems>::Coord
		ExKernelT<ExItems>::calc_centroid(const FacetHandle& _fh){

			HalfedgeHandle& hh = halfedge_handle(_fh);
			HalfedgeHandle& n_hh = next_halfedge_handle(hh);
			HalfedgeHandle& pre_hh = prev_halfedge_handle(hh);

			VertexHandle& vh = vertex_handle(hh);
			VertexHandle& n_vh = vertex_handle(n_hh);
			VertexHandle& pre_vh = vertex_handle(pre_hh);

			return Coord(coord(vh) +
				coord(n_vh)+
				coord(pre_vh) )/3.0;

	}
	
	////////////////////////////////////////////////
	template<class ExItems>
	double
		ExKernelT<ExItems>::getSigmaC()
	{
			double sigmaC = 0.0;
			int num = 0;
			FacetIterator fi = facet_begin();
			for (fi; fi != facet_end(); fi++)
			{
				HalfedgeHandle& hh = fi->halfedge_handle_;
				FacetHandle&    fh = facet_handle(hh);
				Coord& ci = calc_centroid(fh);

				FacetHandles    NeighborFacets;
				getNeighborFaceRing(fh, 2, NeighborFacets);
				for (int j = 0; j < NeighborFacets.size(); j++)
				{
					FacetHandle& Nfh = NeighborFacets[j];
					Coord& cj = calc_centroid(Nfh);
					sigmaC += (ci - cj).length();
					num++;
				}
			}
			sigmaC = sigmaC / num;
			return sigmaC;
		}

	////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::bilateralNormalFilter(std::vector<Normal>& bilateralGuidanceNormal, double param)
	{
			double SigmaC = getSigmaC();
			double MaxSigmaS = param;

			std::vector<double> facetareas;
			FacetIterator fi(facet_begin());
			for (; fi < facet_end(); fi++)
			{
				HalfedgeHandle& hh = fi->halfedge_handle_;
				FacetHandle&    fh = facet_handle(hh);
				facetareas.push_back(calc_facet_area(fh));
			}

			int num = 10;

			for (int i = 0; i < num; i++)
			{
				std::vector<Normal> new_facet_normals;

				for (fi = facet_begin(); fi < facet_end(); fi++)
				{
					HalfedgeHandle& hh = fi->halfedge_handle_;
					FacetHandle&    fh = facet_handle(hh);
					FacetHandles    NeighorFacets;
					getNeighborFaceRing(fh, 2, NeighorFacets);

					Coord& center_fh = calc_centroid(fh);
					Normal& nf = normal(fh);
					Normal  trans(0, 0, 0);
					for (int j = 0; j<NeighorFacets.size(); j++)
					{
						FacetHandle& Nfh = NeighorFacets[j];
						Coord& Ncenter = calc_centroid(Nfh);
						Normal& NorN = normal(Nfh);

						double Wc, Ws;
						Wc = exp(-(center_fh - Ncenter).sqNorm() / (2 * SigmaC*SigmaC));
						Ws = exp(-(nf - NorN).sqNorm() / (2 * MaxSigmaS*MaxSigmaS));

						trans += NorN * Wc*Ws*facetareas[Nfh.idx()];
					}

					if (trans.length() > 1E-7)
					{
						new_facet_normals.push_back(trans.normalize());
					}
					else
					{
						new_facet_normals.push_back(trans);
					}
				}

				for (fi = facet_begin(); fi<facet_end(); fi++)
				{
					FacetHandle fhh = facet_handle(fi->halfedge_handle_);
					fi->normal_ = new_facet_normals[fhh.idx()];
				}
			}

			//-------calculate vertex normal-------
			VertexIterator vit(vertex_begin());
			for (; vit < vertex_end(); vit++)
			{
				Normal n = calc_normal(vertex_handle(vit->halfedge_handle_));
				bilateralGuidanceNormal.push_back(n);
			}
		}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::getNeighborFaceRing(FacetHandle& _fh, int _ring, std::vector<FacetHandle>& _fhs)
	{
			_fhs.push_back(_fh);
			int iteration = 0;
			int verNewNum = _fhs.size();
			int verOldNum = verNewNum - 1;
			int verOldNum1 = verOldNum;

			do{
				verOldNum = _fhs.size();
				for (int i = verOldNum1; i<verNewNum; i++){
					FacetHandle& fh = _fhs[i];
					HalfedgeHandle& hh = halfedge_handle(fh);
					HalfedgeHandle& pre_hh = prev_halfedge_handle(hh);
					HalfedgeHandle& nex_hh = next_halfedge_handle(hh);

					VertexHandle  vhs[3];
					vhs[0] = vertex_handle(hh);
					vhs[1] = vertex_handle(pre_hh);
					vhs[2] = vertex_handle(nex_hh);

					for (int i = 0; i<3; i++){

						HalfedgeHandle& hhv = halfedge_handle(vhs[i]);
						HalfedgeHandle cursor(hhv);

						do{

							FacetHandle fh = facet_handle(cursor);
							if (fh.is_valid() && fh != _fh){
								int j = 0;
								for (; j< _fhs.size(); j++){

									if (fh.idx() == _fhs[j].idx()) break;
								}//end for

								if (j >= _fhs.size()) _fhs.push_back(fh);
							}

							cursor = cw_rotated(cursor);

						} while (hhv != cursor);//end for do while
					}//end for
				}
				verNewNum = _fhs.size();
				verOldNum1 = verOldNum;
				iteration++;
			} while (iteration < _ring);
		}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////_initialization//////////////////////////////////////////
	template<class ExItems>
	bool
		ExKernelT<ExItems>::_initialization(){

			patchNum = facet_number();
			vertexNum = vertex_number();

			if ((patchNum == 0) || (vertexNum == 0)) return false;

			matA = new CSpMatrix(patchNum*FACENUM, vertexNum);
			matATA = new CSpMatrix(vertexNum, vertexNum);

			vectorRz = new double[patchNum*FACENUM];
			vectorXz = new double[vertexNum];
			vectorBz = new double[vertexNum];

			vectorRy = new double[patchNum*FACENUM];
			vectorXy = new double[vertexNum];
			vectorBy = new double[vertexNum];

			vectorRx = new double[patchNum*FACENUM];
			vectorXx = new double[vertexNum];
			vectorBx = new double[vertexNum];

			return true;
		}
	//////////////////////////////void _fillMatrixA//////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::_fillMatrixA(){
			int patchId = 0;


			FacetIterator fi(facet_begin());
			for (; fi<facet_end(); fi++, patchId++){


				HalfedgeHandle& hh = fi->halfedge_handle_;
				FacetHandle&    fh = facet_handle(hh);
				VertexHandle&   vh = vertex_handle(hh);

				int vertexIdArray[FACENUM];
				for (int i = 0; i<FACENUM; i++) {
					vertexIdArray[i] = vh.idx();

					hh = next_halfedge_handle(hh);
					vh = vertex_handle(hh);
				}

				//double coeff1 = (-1.0/((double)FACENUM));
				//double coeff2 = 1.0+coeff1;

				for (int i = 0; i<FACENUM; i++){
					if (vertexIdArray[i]<0) continue;
					for (int j = 0; j<FACENUM; j++){
						if (i == j){
							matA->Set(patchId*FACENUM + j, vertexIdArray[i], 1);//coeff2);
							//matA->Get(patchId*FACENUM+j, vertexIdArray[i], val);
							//std::cout<<val<<" ";
						}
						/*else{
						matA->Set(patchId*FACENUM+j, vertexIdArray[i],coeff1);
						}*/

					}
				}
			}
			CSpMatrix::MulATA(*matA, *matATA);
		}
	/////////////////////////////_factorizeMatrixA///////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::_factorizeMatrixA(){
			if (!LinearSolver::FactorA(*matATA, factorized)) printf("Cannot factorize Matrix!\n");

		}
	/////////////////////////projection///////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::_projection(){//////////////////////projected face
			int number = facet_number();

			int FaceId = 0;
			int n = 0;

			FacetIterator fi(facet_begin());
			for (; fi<facet_end(); fi++, FaceId++){


				HalfedgeHandle& hh = fi->halfedge_handle_;
				FacetHandle&    fh = facet_handle(hh);
				VertexHandle&   vh = vertex_handle(hh);
				Coord& ord = vertex_ref(vh).coord_;

				Normal& norn = normal(fh);
				Coord& before_center = calc_centroid(fh);


				double pp[3], center[3], Overshoot_valueX[3], Overshoot_valueY[3], Overshoot_valueZ[3], p[3];///////////A surface mapping point pp[i]
				///central point center[i] for a surface
				double len[3];
				///Temporary point overshoot_value[i]
				for (int i = 0; i<3; i++) {
					center[i] = 0.0;
					len[i] = 0.0;
				}

				//std::vector<Coord> coordinate,original_coordinate; 

				for (int i = 0; i<FACENUM; i++){
					Coord& varry = vertex_ref(vh).coord_;


					pp[0] = varry.x;
					pp[1] = varry.y;
					pp[2] = varry.z;


					p[0] = (-(pp[2] - before_center[2])*norn.z*norn.x - (pp[1] - before_center[1])*norn.y*norn.x +
						pp[0] * (norn.z*norn.z + norn.y*norn.y) + before_center[0] * norn.x*norn.x);

					p[1] = (-(pp[0] - before_center[0])*norn.x*norn.y - (pp[2] - before_center[2])*norn.z*norn.y +
						pp[1] * (norn.x*norn.x + norn.z*norn.z) + before_center[1] * norn.y*norn.y);

					p[2] = (-(pp[0] - before_center[0])*norn.x*norn.z - (pp[1] - before_center[1])*norn.y*norn.z +
						pp[2] * (norn.x*norn.x + norn.y*norn.y) + before_center[2] * norn.z*norn.z);

					//coordinate.push_back(Coord(p[0],p[1],p[2]));
					//original_coordinate.push_back(varry);
					Coord coordinate = Coord(p[0], p[1], p[2]);
					len[i] = (varry - coordinate).tLength();

					Overshoot_valueX[i] = p[0]; Overshoot_valueY[i] = p[1]; Overshoot_valueZ[i] = p[2];
					center[0] += p[0]; center[1] += p[1]; center[2] += p[2];

					hh = next_halfedge_handle(hh);
					vh = vertex_handle(hh);
				}

				center[0] /= (float)FACENUM;
				center[1] /= (float)FACENUM;
				center[2] /= (float)FACENUM;


				for (int i = 0; i<FACENUM; i++){
					vectorRx[FaceId*FACENUM + i] = Overshoot_valueX[i];//-center[0];
					vectorRy[FaceId*FACENUM + i] = Overshoot_valueY[i];//-center[1];
					vectorRz[FaceId*FACENUM + i] = Overshoot_valueZ[i];//-center[2];
				}

			}
			CSpMatrix::MulATB(*matA, vectorRz, vectorBz);
			CSpMatrix::MulATB(*matA, vectorRy, vectorBy);
			CSpMatrix::MulATB(*matA, vectorRx, vectorBx);
		}
	///////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////_backSubstitution//////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::_backSubstitution(){
			if (!LinearSolver::SolveX(factorized, vectorXz, vectorBz)) printf("Cannot solve!\n");
			if (!LinearSolver::SolveX(factorized, vectorXy, vectorBy)) printf("Cannot solve!\n");
			if (!LinearSolver::SolveX(factorized, vectorXx, vectorBx)) printf("Cannot solve!\n");
		}

	/////////////////////////////_updateResult////////////////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::_updateResult(){

			FacetIterator fi(facet_begin());
			for (; fi<facet_end(); fi++){
				HalfedgeHandle& hh = fi->halfedge_handle_;
				FacetHandle&    fh = facet_handle(hh);
				VertexHandle&   vh = vertex_handle(hh);

				for (int i = 0; i<FACENUM; i++){

					vertex_ref(vh).coord_.x = vectorXx[vh.idx()];
					vertex_ref(vh).coord_.y = vectorXy[vh.idx()];
					vertex_ref(vh).coord_.z = vectorXz[vh.idx()];

					hh = next_halfedge_handle(hh);
					vh = vertex_handle(hh);

				}

			}

		}
	///////////////////////////////////////////////////////////////
	//////////////global vertex updating//////////////////////
	/////////////////////////////DoComputing ///////////////////////////////////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::DoComputing(int iVertexItr, float fMaxSigma, float fMaxSigmaC, int iFacetItr){
			///////////DoComputing Rolling Guidance Filter and Surface-from-Gradients

			adjustFaceNormal_RGF(fMaxSigma, fMaxSigmaC, iFacetItr);
			globalVertexUpdating(iVertexItr);
		}
	///////////////////////////////////////////////////////////////
	//////////////global vertex updating//////////////////////
	template<class ExItems>
	void
		ExKernelT<ExItems>::globalVertexUpdating(int iVertexItr){
			_initialization();
			_fillMatrixA();
			_factorizeMatrixA();
			for (int i = 0; i<iVertexItr; i++){
				_projection();
				_backSubstitution();
				_updateResult();
			}

			if (factorized){
				LinearSolver::DeleteF(factorized);
				factorized = NULL;
			}
		}

} /// namespace



