///////////////////////////////////////////////////////////////////////////////
#ifndef __WMQ_MESH_KERNELT_H__
#define __WMQ_MESH_KERNELT_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <assert.h> 
#include "Vector3T.h" 

#include "..\LinearSolver\linearsolver.h"

namespace MeshN {


	///////////////////////////////////////////////////////////////////////////////
	// KernelT: Mesh kernel, which includes basic operators for low-data managing
	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	class KernelT {
	public: // definitions of types 
		// definition of this class
		typedef typename KernelT<Items>         This;

		// definitions of non-topologic types
		typedef typename Items::Scalar          Scalar;//数量、标量
		typedef typename Items::Coord           Coord;//坐标
		typedef typename Items::Normal          Normal;//向量
		typedef typename Items::Color           Color;//颜色
		typedef typename Items::TexCoord        TexCoord;//纹理坐标

		// definitions of topologic items
		typedef typename Items::Halfedge        Halfedge;//半边
		typedef typename Items::Vertex          Vertex;//顶点
		typedef typename Items::Edge            Edge;//整边
		typedef typename Items::Facet           Facet;//面

		// definitions of various handles
		typedef typename Items::HalfedgeHandle  HalfedgeHandle;//半边的句柄，相当于指针
		typedef typename Items::VertexHandle    VertexHandle;//顶点的句柄
		typedef typename Items::EdgeHandle      EdgeHandle;//整边的句柄
		typedef typename Items::FacetHandle     FacetHandle;//面的句柄

		// handle arrays
		typedef typename std::vector<VertexHandle>   VertexHandles;
		typedef typename std::vector<EdgeHandle>     EdgeHandles;
		typedef typename std::vector<FacetHandle>    FacetHandles;

		// iterator types
		typedef typename std::vector<Vertex>::iterator  VertexIterator;
		typedef typename std::vector<Edge>::iterator    EdgeIterator;
		typedef typename std::vector<Facet>::iterator   FacetIterator;
		// constant iterator types
		typedef typename std::vector<Vertex>::const_iterator  ConstVertexIterator;
		typedef typename typename std::vector<Edge>::const_iterator    ConstEdgeIterator;
		typedef typename std::vector<Facet>::const_iterator   ConstFacetIterator; 


	public: // convenient functions  
		// definitions of following funcs in "ConvenientFuncs.cpp" 
		// half edge handle <--> ?? handles 
		inline HalfedgeHandle  next_halfedge_handle(const HalfedgeHandle& _hh);
		inline HalfedgeHandle  prev_halfedge_handle(const HalfedgeHandle& _hh);
		inline VertexHandle    vertex_handle(const HalfedgeHandle& _hh);
		inline EdgeHandle      edge_handle(const HalfedgeHandle& _hh);
		inline FacetHandle     facet_handle (const HalfedgeHandle& _hh);
		inline HalfedgeHandle  opposite_halfedge_handle(const HalfedgeHandle& _hh);

		// ?? handles <--> halfedge handle: get halfedge-handle from v/e/f
		inline HalfedgeHandle  halfedge_handle(const VertexHandle& _vh);
		inline HalfedgeHandle  halfedge_handle(const EdgeHandle& _eh, int _i);
		inline HalfedgeHandle  halfedge_handle(const FacetHandle& _fh);

		inline HalfedgeHandle ccw_rotated(const HalfedgeHandle& _hh); 
		inline HalfedgeHandle  cw_rotated(const HalfedgeHandle& _hh);

		// handles <--> refereces
		inline Halfedge&  halfedge_ref(const HalfedgeHandle& _hh); ////&
		inline Vertex&    vertex_ref(const VertexHandle& _vh);  
		inline Edge&      edge_ref(const EdgeHandle& _eh); 
		inline Facet&     facet_ref(const FacetHandle& _fh); 

	public: // Manage status  
		// definitions of following funcs in "StatusManager.cpp" 
		inline void set_deleted(const VertexHandle& _vh, bool _b);
		inline void set_deleted(const EdgeHandle& _eh, bool _b);
		inline void set_deleted(const FacetHandle& _fh, bool _b);
		inline void set_deleted(const HalfedgeHandle& _hh, bool _b);

		inline bool is_deleted(const VertexHandle& _vh);
		inline bool is_deleted(const HalfedgeHandle& _hh);
		inline bool is_deleted(const EdgeHandle& _eh);
		inline bool is_deleted(const FacetHandle& _fh);

		// and so on
		// ........



	public: // ctor and dtor and assign 
		// definitions of following funcs in "StaticKernelT.cpp" 
		// ctor and dtor
		//KernelT() {}

		~KernelT() { clear(); }
		KernelT() : matA(NULL)
			, matATA(NULL)
			, vectorRz(NULL)
			, vectorBz(NULL)
			, vectorXz(NULL)
			, factorized(NULL)
			, vectorRy(NULL)
			, vectorRx(NULL)
			, vectorXy(NULL)
			, vectorBy(NULL)
			, vectorXx(NULL)
			, vectorBx(NULL)
		{}
		// reset this class
		inline void clear(void);
		// obj assigning 
		inline This& operator=(const This& _rhs); 

	public: // building data structure 
		// definitions of following funcs in "BuildKernel.cpp" 
		// to add a vertex into mesh
		inline VertexHandle add_vertex(void);
		inline VertexHandle add_vertex(const Coord& _coord);

		// to add a facet into mesh
		inline FacetHandle add_facet(const std::vector<VertexHandle>& _vhs); 
		inline FacetHandle add_facet(const VertexHandle& _vh0, 
			const VertexHandle& _vh1,
			const VertexHandle& _vh2 ); 

		// to add edge to mesh
		inline EdgeHandle add_edge(const VertexHandle& _vh0,
			const VertexHandle& _vh1);
	private: 
		// find halfedge <_vh0, _vh1>, if find no, return invalid_handle
		inline HalfedgeHandle _find_halfedge_handle(const VertexHandle& _vh0,
			const VertexHandle& _vh1);

		// correct order error when a facet is added to mesh
		inline void _fix_halfedge_order(const HalfedgeHandle& _hh0,
			const HalfedgeHandle& _hh1);
		// adhere a(only) boundary_halfedge to given vertex 
		inline bool _adjust_to_boundary_halfedge(const VertexHandle& _vh);


	public: // checkings 
		// definitions of following funcs in "StaticKernelT.cpp" 
		// boundary? 
		inline bool is_boundary_halfedge(const HalfedgeHandle& _heh);
		inline bool is_boundary_edge(const EdgeHandle& _eh);
		inline bool is_boundary_vertex(const VertexHandle& _vh);
		inline bool is_boundary_facet(const FacetHandle& _fh);

		// number(exclude deleted ones)
		inline int vertex_number(void); 
		inline int edge_number(void);
		inline int facet_number(void);

		// total size
		inline int vertex_size(void) { return vertices_.size(); }
		inline int halfedge_size(void) { return edges_.size()<<1; }
		inline int edge_size(void) { return edges_.size(); }
		inline int facet_size(void) { return facets_.size(); }

		// coordinates of vertices
		inline Coord& coord(const VertexHandle& _vh);

	public:
		// iterators
		inline VertexIterator vertex_begin(void) { return vertices_.begin(); }
		inline VertexIterator vertex_end(void)   { return vertices_.end(); }

		inline EdgeIterator   edge_begin(void)   { return edges_.begin(); }
		inline EdgeIterator   edge_end(void)     { return edges_.end(); }

		inline FacetIterator  facet_begin(void)  { return facets_.begin(); }
		inline FacetIterator  facet_end(void)    { return facets_.end(); }


	public:
		inline Scalar get_average_edge_length()const{//取平均边长

			return averaged_edge_length_;
		}

		inline void getMeshBox(Coord& _min, Coord& _max)const{

			_min = min_;
			_max = max_;
		}

	public:
		inline void set_average_edge_length(Scalar _edgelength){averaged_edge_length_ = _edgelength;}//设置平均边长
		inline void setMeshBox(const Coord& _min, const Coord& _max){//设置网格包围球
			max_ = _max;
			min_ = _min;
		}

		inline bool calcMeshBox();//计算网格包围球

	public: // Data body  
		std::vector<Vertex>  vertices_;
		std::vector<Edge>    edges_;
		std::vector<Facet>   facets_;

	private:
		Scalar averaged_edge_length_;//平均边长
		Coord min_, max_;///包围球最小最大坐标

	

	public:
		std::vector<Vertex> original_vertices_;
		std::vector<Edge> original_edges_;
		std::vector<Facet> original_facets_;
	 public:
		 CSpMatrix* matA, *matATA;//input matrix A which is M X N (M>N) 
		 // and its least square ATA which is N X N 
		 double* vectorRz, *vectorBz, *vectorXz;        //vectorR: right hand side vector (M X 1)
		 // vectorB: right hand side vector which equals to AT X vectorR (N X 1);
		 // vectorX: unknown depth for all vertex which is N X 1
		 void* factorized;// factorized matrix (for direct solver)

		 int patchNum, vertexNum; // total patch number and vertex number

		 double* vectorRy, *vectorXy, *vectorBy;
		 double* vectorRx, *vectorXx, *vectorBx;
	}; 


} /// namespace

#endif 

