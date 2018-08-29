#include "BaseHandle.h"
//////////////implementation/////////////////
namespace MeshN{

bool MeshN::BaseHandle::operator== (const MeshN::BaseHandle& _rhs) const { 
	return idx_ == _rhs.idx_; 
}

///////////////////////////////////////////////////////////////////////////////
bool MeshN::BaseHandle::operator!= (const MeshN::BaseHandle& _rhs) const { 
	return idx_ != _rhs.idx_; 
}

///////////////////////////////////////////////////////////////////////////////
bool MeshN::BaseHandle::operator<  (const MeshN::BaseHandle& _rhs) const { 
	return idx_ < _rhs.idx_; 
}

///////////////////////////////////////////////////////////////////////////////
bool MeshN::BaseHandle::operator>  (const MeshN::BaseHandle& _rhs) const { 
	return idx_ > _rhs.idx_; 
}

}