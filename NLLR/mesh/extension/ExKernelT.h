
#ifndef __WMQ_MESH_EXTENSION_KERNEL_H__
#define __WMQ_MESH_EXTENSION_KERNEL_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

///////////////////////////////////////////////////////////////////////////////
#include "../core/KernelT.cpp" 
#include "time.h"
#include <string>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include <linearsolver.h>

using namespace Eigen;
using namespace std;
#define   weimq_e 2.71828
#define   FACENUM 3
namespace MeshN { 

	template <class ExItems>
	class ExKernelT : public MeshN::KernelT<ExItems> {
	public:  
		typedef typename MeshN::KernelT<ExItems>       Kernel; 
		typedef typename ExKernelT<ExItems>            This; 

		typedef typename Kernel::Scalar            Scalar;
		typedef typename Kernel::Coord             Coord;
		typedef typename Kernel::Normal            Normal;  
		typedef typename Kernel::Color             Color;
		typedef typename Kernel::TexCoord          TexCoord;

		typedef typename Kernel::Halfedge          Halfedge;
		typedef typename Kernel::Vertex            Vertex;
		typedef typename Kernel::Edge              Edge;
		typedef typename Kernel::Facet             Facet;

		typedef typename Kernel::HalfedgeHandle    HalfedgeHandle;
		typedef typename Kernel::VertexHandle      VertexHandle;
		typedef typename Kernel::EdgeHandle        EdgeHandle;
		typedef typename Kernel::FacetHandle       FacetHandle; 

		typedef typename Kernel::VertexHandles     VertexHandles;
		typedef typename Kernel::EdgeHandles       EdgeHandles;
		typedef typename Kernel::FacetHandles      FacetHandles;

		typedef typename Kernel::VertexIterator       VertexIterator;
		typedef typename Kernel::EdgeIterator         EdgeIterator;
		typedef typename Kernel::FacetIterator        FacetIterator; 
		typedef typename Kernel::ConstVertexIterator  ConstVertexIterator;
		typedef typename Kernel::ConstEdgeIterator    ConstEdgeIterator;
		typedef typename Kernel::ConstFacetIterator   ConstFacetIterator; 

	public:
		ExKernelT();
		~ExKernelT();


	public: 
		// normals of facets
		inline Normal normal(const FacetHandle& _fh);
		inline Normal calc_normal(const FacetHandle& _fh);
		void   update_facet_normals(void);
		// normals of vertices
		inline Normal normal(const VertexHandle& _vh);
		inline Normal calc_normal(const VertexHandle& _vh);
		inline Normal calc_normal_max(const VertexHandle& _vh);
		void   update_vertex_normals(void);//
		void   update_vertex_normals_max(void);
		Scalar calc_facet_area(const FacetHandle& _fh);
		void   update_area();
		Scalar get_area(const FacetHandle& _fh);
		double calc_edge_length(const EdgeHandle& _eh);

		void   update_normals(void);

		//edge length...
		void update_edge_length(void);
		void output_to_file();
		void output_to_file(char* filename);
		void getNeighborRing(VertexHandle& _vh, int _ring, std::vector<VertexHandle>& _vhs);//取顶点的ring 环的所有点
		void getNeighborVertices(VertexHandle& _vh, int _verticesNum, std::vector<VertexHandle>& NeighborVertices);
		void getNeighborVertices_withRingIndex(VertexHandle& _vh, int _verticesNum, std::vector<VertexHandle>& NeighborVertices, std::vector<int>& ringIdx);
		void getNeighborFaceN1(FacetHandle& _fh, std::vector<FacetHandle>& _fhs);//getting faces sharing edges with _fh
		void getNeighborFaceN2(FacetHandle& _fh, std::vector<FacetHandle>& _fhs);//getting faces sharing common vertices with _fh

		Coord calc_centroid(const FacetHandle& _fh);

		void getNeighborFaceRing(FacetHandle& _fh, int _ring, std::vector<FacetHandle>& _fhs);

		//////////////////////////////////////////////////////////////////////////
		double computeSig(double nsig);
		void updateNormalIterativeRegularization(double w);

		//////////////////////////////////////////////////////////////////////////
		void meshInit();
		void findSimilarPatch3(int _candidateSize, int _t); //normal region covariance + guidance
		void patchNormalAlignment(std::vector<Normal> guidanceNormals, int K);
		void getRegionNormalCovarianceMat_new(std::vector<Normal> guidanceNormals, int K);
		void low_rank_recover(int _t, double c1, double nsig, int K);
		void low_rank_ssc_C(MatrixXd Y, double c1, double nsig, MatrixXd &output, MatrixXd &W);
		void low_rank_ssc_C_TG(MatrixXd Y, double c1, double nsig, MatrixXd &output, MatrixXd &W);
		MatrixXd admmAXB_TG(MatrixXd A, MatrixXd B, MatrixXd X, double rho, MatrixXd thr, int R);
		void DC(MatrixXd D, double rho, MatrixXd thr, MatrixXd T0, double a, MatrixXd &X, MatrixXd &T);
		void compute_normals_NLLR();
		void adjustVertexCoord(int iterations);
		void bilateralNormalFilter(std::vector<Normal>& bilateralGuidanceNormal, double param); //Based on bilateral normal filtering (vertex normal)
		void calcUpdatedVertexNormal(std::vector<Normal>& updatedVertexNormals);
		double getSigmaC();

		/////////global vertex updating//////////////////////////////////////////
		void globalVertexUpdating(int iVertexItr);
		bool _initialization();
		void _fillMatrixA();
		void _factorizeMatrixA();
		void _projection();
		void _backSubstitution();
		void _updateResult();
		void DoComputing(int iVertexItr, float fMaxSigma, float fMaxSigmaC, int iFacetItr);
		//////////////////////////////////////////////////////////////////////////
		void update_flag(void);
		void get_original_vertex();
		void get_original_facet();

		//////////////////////////// global vertex updating 2 //////////////////////////////////
		//get point handles for face
		void getVextexHandlesForFace(FacetHandle &_fh, std::vector<VertexHandle> &_vhs);
		//get point coordinates for face
		void getPointsForFace(FacetHandle &_fh, std::vector<Coord> &_vhs);
		//get all the faces that contain the point
		void getFacesForPoint(VertexHandle vh, std::vector<FacetHandle> &faceHVector);
		//calculate the normals of each face and store it in an array
		void getFaceNormal(std::vector<Normal> &normalVector);
		//calculate the area of each face and store it in an array
		void getFaceAera(std::vector<double> &areaVector);

		//update the vertex position with the quadratic optimization
		void UpdateVertexPosition(std::vector<Normal> &original_facet_normals, std::vector<Normal> &filtered_normals);
		//void getNeighborFaceRing(FacetHandle& _fh, int _ring, std::vector<FacetHandle>& _fhs); %lzhu

		//calculate the new gradient vectors with the filtered normals
		void calculateNewGradientVector(std::vector<Normal> &original_normals, std::vector<Normal> &filtered_normals, std::vector<std::vector<Coord> > &rotated_gradients);
		//get the 3x3 rotation matrix for each triangle face
		void getRotationMatrixForEachFace(Normal &original_normal, Normal &filtered_normal, Eigen::Matrix3d &rotation_matrix);
		//calculate the divergence of each vertex
		void calculateDivergence(std::vector<std::vector<Coord> >&rotated_gradients, std::vector<std::vector<double> > &divergenceVector);

		void getPoissonCoefficientMatrix(Eigen::SparseMatrix<double> &p_Coefficient_Matrix);

		//get the rotation matrix in the Esmooth term
		void calculateFaceRotationMatrix(Coord rotationAxisVec, Normal firstFaceNormal,
			Normal secondFaceNormal, Eigen::Matrix3d &rotationMatrix);
		//update coefficient matrix of the Esmooth item with quadratic terms
		void updateCoefficientMatrixWithQuadraticTerm(int faceHandleIndex, double edgeLength, double weightedValue,
			Eigen::SparseMatrix<double> &s_Coefficient_Matrix);
		//update coefficient matrix of the Esmooth item with cross terms
		void updateCoefficientMatrixWithCrossTerm(int firstFaceHandleIndex, int secondFaceHandleIndex,
			int conjunctEdgeIndex, double weightedValue, Eigen::SparseMatrix<double> &s_Coefficient_Matrix);
		//get the coefficient matrix of the Esmooth term
		void getSmoothCoefficientMatrix(Eigen::SparseMatrix<double> &s_Coefficient_Matrix);

		//solve the sparse linear equation
		void solveLinearEquation(Eigen::SparseMatrix<double> &t_Coefficient_Matrix,
			std::vector<std::vector<double> > &divergenceVector, Eigen::MatrixXd &coordinatesVector);

		void solveLinearEquation_CG(Eigen::SparseMatrix<double> &t_Coefficient_Matrix,
			std::vector<std::vector<double> > &divergenceVector, Eigen::MatrixXd &coordinatesVector);
		//////////////////////////// global vertex updating 2 //////////////////////////////////

	public:
		bool has_Normal(){return isNormal_;}
		void set_isNormal(bool tf){isNormal_ = tf;}
		bool has_Area(){return isArea_;}
		void set_isArea(bool tf){isArea_ = tf;}

	private:
		bool isNormal_;
		bool isArea_;
	};

} /// namespace

#endif



