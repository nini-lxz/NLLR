/////////////////////////////////////////////////////////////////////////////// 
#ifndef __WMQ_MESH_EXTENSION_ITEMS_H__
#define __WMQ_MESH_EXTENSION_ITEMS_H__ 

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "../core/InternalItems.h"
namespace MeshN {


	///////////////////////////////////////////////////////////////////////////////
	class ExItems : public MeshN::InternalItems{
	public:   
		struct Vertex : public MeshN::InternalItems::Vertex {  
			Normal   normal_; 
			Normal         oriNormal_;
		};


		///===================================================
		struct Edge : public MeshN::InternalItems::Edge {
			Scalar   length_;
		};


		///===================================================
		struct Facet : public MeshN::InternalItems::Facet {

			Normal         normal_; 
			Normal         oriNormal_;
			Coord          centroid_;
			Scalar         area_;
			//add others

		}; 

	};


} //namespace 


#endif

