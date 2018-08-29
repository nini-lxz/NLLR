
#ifndef __WMQ_MESH_INTERNAL_HANDLES_H__
#define __WMQ_MESH_INTERNAL_HANDLES_H__ 

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "BaseHandle.cpp"
namespace MeshN {
	struct InternalHandles {
		// Handle for a vertex entity
		struct VertexHandle : public BaseHandle {
			explicit VertexHandle(int _idx=-1) : BaseHandle(_idx) {}
		};

		// Handle for a halfedge entity
		struct HalfedgeHandle : public BaseHandle {
			explicit HalfedgeHandle(int _idx=-1) : BaseHandle(_idx) {}
		};

		// Handle for an edge entity
		struct EdgeHandle : public BaseHandle {
			explicit EdgeHandle(int _idx=-1) : BaseHandle(_idx) {}
		};

		// Handle for a face entity
		struct FacetHandle : public BaseHandle {
			explicit FacetHandle(int _idx=-1) : BaseHandle(_idx) {}
		};  
	}; // InternalHandles

}  // namespace


#endif