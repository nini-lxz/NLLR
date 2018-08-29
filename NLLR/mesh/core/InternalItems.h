/////////////////////////////////////////////////////////////////////////////// 
#ifndef __WMQ_MESH_INTERNAL_ITEMS_H__
#define __WMQ_MESH_INTERNAL_ITEMS_H__ 

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

///////////////////////////////////////////////////////////////////////////////
#include "Vector3T.h"
#include "InternalHandles.h"
#include "InternalStatus.h"
#include <stdio.h>
#include <iostream>
#define NN 20
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
using namespace Eigen;

namespace MeshN {

	class InternalItems {
	public:
		// definitions of types
		typedef float           Scalar; 
		typedef MathN::Vec3f    Coord; 
		typedef MathN::Vec3f    Normal;
		typedef MathN::Vec3f    Color; 
		typedef MathN::Vec3f    TexCoord;

		// handle types
		typedef InternalHandles::HalfedgeHandle  HalfedgeHandle;
		typedef InternalHandles::VertexHandle    VertexHandle;
		typedef InternalHandles::EdgeHandle      EdgeHandle;
		typedef InternalHandles::FacetHandle     FacetHandle;

	public:
		// topologic --- halfedge ////////////////////
		struct Halfedge {
			HalfedgeHandle    next_halfedge_handle_;
			HalfedgeHandle    prev_halfedge_handle_;
			VertexHandle      vertex_handle_;
			FacetHandle       facet_handle_;
			Status            status_; 
		};

		// topologic --- vertex //////////////////////
		struct Vertex { 
			HalfedgeHandle    halfedge_handle_; 
			Status            status_;
			Coord             oriCoord_; //original vertex position
			Coord             coord_;
			//int               similarIdx_[NN];
			MatrixXd          coMatrix_;
			MatrixXd          guiCoMatrix_;
			MatrixXd          rotMat_;
			MatrixXd          patch_aligned_;
			MatrixXd          patch_aligned_gui_;
			std::vector<int>  similarIdx_;
			double            featureVal_;
			double            guidanceVal_;
			MathN::Vec3f      sumNormal_;
			MathN::Vec3f      weights_;
			std::vector<VertexHandle> neighborVertices_; //xz
		};

		// topologic --- edge ////////////////////////
		struct Edge {
			Halfedge          halfedges_[2];//two halfedges
			Status            status_;  
		};

		// topologic --- facet ///////////////////////
		struct Facet {
			HalfedgeHandle    halfedge_handle_;
			Status            status_;  
		}; 

	}; // InternalItems


}  // namespace 


#endif

