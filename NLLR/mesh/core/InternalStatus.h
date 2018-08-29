
#ifndef __WMQ_MESH_INTERNAL_STATUS_H__
#define __WMQ_MESH_INTERNAL_STATUS_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

namespace MeshN { 

	enum StatusBits { 
		NONE      =   0, // none
		DELETED   =   1, // data is deleted
		LOCKED    =   2, // data is locked
		HIDDEN    =   4, // data is hidden
		MODIFIED  =   8, // data is modified
		SELECTED  =  16, // data is selected.
		FEATURE   =  32, // data is or belongs to a feature.
		TAGGED    =  64, // data is tagged
		RESERVED  = 128  // this bit is reserved for other use
	};


	/////////////////////////////////////////////////////////////////////////////// 
	class Status {
	public:
		// status data body
		typedef unsigned int StatusByte;
		typedef unsigned int StatusBit;

	public:
		// ctor--> clear status-bits
		Status() : status_byte_(0) {}
		// reset status bits to init (0)
		inline void reset(void) { status_byte_ = 0; }

	public:
		// is a certain bit set ?
		inline bool is_bit_set(StatusBit _s) const { return (status_byte_&_s)>0; }
		// set a certain bit
		inline void set_bit(StatusBit _s) { status_byte_ |= _s; }
		// unset a certain bit
		inline void unset_bit(StatusBit _s) { status_byte_ &= ~_s; }
		// set or unset a certain bit
		inline void change_bit(StatusBit _s, bool _b) {  
			(_b) ? (status_byte_ |= _s) : (status_byte_ &= ~_s);  }

		// get whole status
		inline StatusBit bits(void) const { return status_byte_; }
		// set whole status at once
		inline void set_bits(StatusBit _bits) { status_byte_ = _bits; }

	public:
		// is deleted ?
		inline bool is_deleted(void) const { return is_bit_set(DELETED); }
		// set deleted
		inline void set_deleted(bool _b) { change_bit(DELETED, _b); }

		// is locked ?
		inline bool is_locked(void) const { return is_bit_set(LOCKED); }
		// set locked
		inline void set_locked(bool _b) { change_bit(LOCKED, _b); }

		// is hidden ?
		inline bool is_hidden(void) const { return is_bit_set(HIDDEN); }
		// set hidden
		inline void set_hidden(bool _b) { change_bit(HIDDEN, _b); }

		// is modified ?
		inline bool is_modified(void) const { return is_bit_set(MODIFIED); }
		// set modified
		inline void set_modified(bool _b) { change_bit(MODIFIED, _b); }

		// is selected ?
		inline bool is_selected(void) const { return is_bit_set(SELECTED); }
		// set selected
		inline void set_selected(bool _b) { change_bit(SELECTED, _b); }

		// is feature ?
		inline bool is_feature(void) const { return is_bit_set(FEATURE); }
		// set feature
		inline void set_feature(bool _b) { change_bit(FEATURE, _b); }

		// is tagged ?
		inline bool is_tagged(void) const { return is_bit_set(TAGGED); }
		// set tagged
		inline void set_tagged(bool _b) { change_bit(TAGGED, _b); }

		/// is reserved ?
		inline bool is_reserved(void) const { return is_bit_set(RESERVED); }
		/// set reserved
		inline void set_reserved(bool _b) { change_bit(RESERVED, _b); }

	private: 
		// status data
		StatusByte  status_byte_;
	}; // class Status


}  // namespace


#endif 
