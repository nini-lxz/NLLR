
#ifndef __WMQ_MESH_BASE_HANDLE_H__
#define __WMQ_MESH_BASE_HANDLE_H__ 

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

namespace MeshN {

	class BaseHandle { 
	private:
		int idx_; 

	public: 
		inline explicit BaseHandle(int _idx = -1) : idx_(_idx) {}//idx_=-1,Constructor

		inline void reset(void) { idx_ = -1; }//

		inline void set(int _idx) { idx_ = _idx; }//
		inline int  idx(void) const { return idx_; }

		inline bool is_valid(void) const { return idx_ != -1; }

		//compare
		inline bool operator==(const BaseHandle& _rhs) const;
		inline bool operator!=(const BaseHandle& _rhs) const;
		inline bool operator< (const BaseHandle& _rhs) const;
		inline bool operator> (const BaseHandle& _rhs) const;

		inline void increment() { ++idx_; }// ×ÔÔöº¯Êı
		inline void decrement() { --idx_; }
	};

}  // namespace

#endif