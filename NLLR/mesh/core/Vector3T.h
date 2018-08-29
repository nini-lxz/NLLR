
/************************************************************************/
/* 
*3d vector template class with a set of basic operators
*by Wei Mingqiang (mingqiang.wei@gmail.com)
*13:56:50.2015.7.16
*all rights reserved                                                                    */
/************************************************************************/

#ifndef __WMQ_MATH_VECTOR3T_H__
#define __WMQ_MATH_VECTOR3T_H__

/* #pragma warning (disable: 4786) */

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>
#include <assert.h>
#include <math.h> 

#ifndef _PI
#define _PI 3.1415926535897932384626433f
#endif

namespace MathN 
{

	///////////////////////////////////////////////////////////////////////////////
	// Declaration of template class for 3D vectors
	///////////////////////////////////////////////////////////////////////////////
	template <typename T>            //T: scalar type of data
	class Vector3T 
	{
	public: 
		typedef T            Scalar;  //type of Scalar;
		typedef Vector3T<T>  Vector;  //type of the class

	public:
		union 
		{
			T data_[3];  //array-form
			union 
			{ //vector-form
				struct { T x, y, z; };  //coordinate: (x,y,z)
				struct { T u, v, w; };  //barycenter coordinate: (u,v,w)
				struct { T r, g, b; };  //RGB color: (r,g,b)
				struct { T A, B, C; };  
			};
		};

	public: 
		inline Vector3T(void);
		inline Vector3T(const T& _s) { x=y=z=_s; }
		inline Vector3T(const T& _x, const T& _y, const T& _z);
		inline Vector3T(const T  _rhs[3]);
		inline Vector3T(const Vector& _rhs);  //copy constructor

		virtual ~Vector3T() {}

		inline Vector& operator=(const Vector& _rhs);   //assignment
		inline Vector& set(const T& _v0, const T& _v1, const T& _v2); 

	public: 
		static inline int dim()      { return 3; }           //get dimension of vector 
		static inline size_t byte()  { return sizeof(T); }   //bytes of T
		static inline size_t bytes() { return 3*sizeof(T); } //bytes of data

		inline T& operator[](int _i);
		inline const T& operator[](int _i) const; 

		// usage: Vector tt; T p = *(tt); OR *(tt+2) = 48;
		inline operator T*(void) { return data_; }
		inline operator const T*(void) const { return data_; }

	public:
		// compare relation between 2 vectors
		inline bool operator==(const Vector& _rhs) const; 
		inline bool operator!=(const Vector& _rhs) const; 
		inline bool operator< (const Vector& _rhs) const; 
		inline bool operator> (const Vector& _rhs) const;
		inline bool operator<<(const Vector& _rhs) const; 
		inline bool operator>>(const Vector& _rhs) const;

		// vector-scalar self- multiply & division
		inline Vector& operator*=(const T& _s);
		inline Vector& operator/=(const T& _s);

		// vector-scalar multiply & division
		inline Vector operator*(const T& _s) const;
		inline Vector operator/(const T& _s) const;

		// vector-vector self- operators
		inline Vector& operator+=(const Vector& _rhs);  // (+=)
		inline Vector& operator-=(const Vector& _rhs);  // (-=) 
		inline Vector& operator|=(const Vector& _rhs);  // (*=) 
		inline Vector& operator/=(const Vector& _rhs);  // (/=)

		// vector-vector operators
		inline Vector operator+(const Vector& _rhs) const;  // (+)
		inline Vector operator-(const Vector& _rhs) const;  // (-)
		inline Vector operator|(const Vector& _rhs) const;  // (|)--vector
		inline Vector operator/(const Vector& _rhs) const;  // (/) 
		inline Vector operator-(void) const;    // unitary minus

		// dot- and cross- products
		inline T    operator*(const Vector& _rhs) const;  // dot-product
		inline Vector operator%(const Vector& _rhs) const;  // cross-product

	public:
		// norm & square norm
		inline T norm(void)     const; // Euclidean norm
		inline T sqNorm(void)   const; // squared Euclidean norm
		inline T length(void)   const { return norm(); }
		inline T sqLength(void) const { return sqNorm(); }
		inline T tLength(void)   const { return norm(); }

		// max, min, mean 
		inline T maxx(void) const;  //max
		inline T minn(void) const;  //min
		inline T mean(void) const;  //mean 

		// minimize & maximize & normalize 
		inline Vector& minimize(const Vector& _rhs);  //minimize 
		inline Vector& maximize(const Vector& _rhs);  //maximize 
		inline Vector& normalize(void);               //normalize

		// vectorize: set components of vector into same value: _s
		inline Vector& vectorize(const T& _s);
		inline Vector& zeros(void) { return vectorize((T)0); }
		inline Vector& ones(void)  { return vectorize((T)1); }

		inline void rotate(Vector& _vec, const T& _angle);
		inline void rotatePoint(Vector& _pnt, const T& _angle);

		inline T angle(const Vector& vec) 
		{
			float c = x*vec.x + y*vec.y + z*vec.z;

			if (length()==0) return 0;
			if (vec.length()==0) return 0;
			c = c/ ( length()*vec.length() );///  
			T a = (float)acos(c);
			return a*57.2958f;  //return a*180/3.14159;
		}

		// perform funtion: T operator()(T) on components of vector
		// usage: eg. "Vector vec; vec.apply(fabs);", where 'fabs' is common math 
		// function in <math.h>, absolutizing the components of vector
		template <typename Func>
		inline Vector& apply(const Func& _func) 
		{
			data_[0] = _func(data_[0]);
			data_[1] = _func(data_[1]);
			data_[2] = _func(data_[2]);

			return *this;
		}

	}; /// Vector3T


	///////////////////////////////////////////////////////////////////////////////
	// Implementations of member functions of Vector3T
	///////////////////////////////////////////////////////////////////////////////

	template <typename T>	//default ctor
	Vector3T<T>::Vector3T(void) : x(0), y(0), z(0) {
		//data_[0] = data_[1] = data_[2] = (T)0;
	}


	template <typename T>
	Vector3T<T>::
		Vector3T(const T& _v0, const T& _v1, const T& _v2) : x(_v0), y(_v1), z(_v2)
	{}


	template <typename T>	                                       //array input
	Vector3T<T>::Vector3T(const T _rhs[3]) {
		memcpy(data_, _rhs, 3*sizeof(T));
	}


	template <typename T>	                                             //clone
	Vector3T<T>::Vector3T(const Vector& _rhs) {
		memcpy(data_, _rhs.data_, 3*sizeof(T));
	}


	template <typename T>                                               //assign
	Vector3T<T>& 
		Vector3T<T>::operator=(const Vector& _rhs) { 
			memcpy(data_, _rhs.data_, 3*sizeof(T)); 
			return *this;
	}


	template <typename T>	                                            //assign
	Vector3T<T>& 
		Vector3T<T>::set(const T& _v0, const T& _v1, const T& _v2) {
			data_[0] = _v0;   data_[1] = _v1;   data_[2] = _v2; 
			return *this;
	}


	template <typename T> 
	T& Vector3T<T>::operator[](int _i) { 
		assert( 0<=_i && _i<3 );
		return data_[_i]; 
	}


	template <typename T> 
	const T& Vector3T<T>::operator[](int _i) const { 
		assert( 0<=_i && _i<3 );
		return data_[_i]; 
	}


	template <typename T> 
	bool Vector3T<T>::operator==(const Vector& _rhs) const {
		if (data_[0] != _rhs.data_[0]) return false;
		if (data_[1] != _rhs.data_[1]) return false;
		if (data_[2] != _rhs.data_[2]) return false;
		return true; 
	}


	template <typename T> 
	bool Vector3T<T>::operator!=(const Vector& _rhs) const {
		if (data_[0] != _rhs.data_[0]) return true;
		if (data_[1] != _rhs.data_[1]) return true;
		if (data_[2] != _rhs.data_[2]) return true;
		return false; 
	}


	template <typename T> 
	bool Vector3T<T>::operator<(const Vector& _rhs) const {
		if (data_[0] != _rhs.data_[0])
			return (data_[0] < _rhs.data_[0]);
		if (data_[1] != _rhs.data_[1])
			return (data_[1] < _rhs.data_[1]);
		if (data_[2] != _rhs.data_[2])
			return (data_[2] < _rhs.data_[2]);
		return false;
	}


	template <typename T> 
	bool Vector3T<T>::operator>(const Vector& _rhs) const {
		if (data_[0] != _rhs.data_[0])
			return (data_[0] > _rhs.data_[0]);
		if (data_[1] != _rhs.data_[1])
			return (data_[1] > _rhs.data_[1]);
		if (data_[2] != _rhs.data_[2])
			return (data_[2] > _rhs.data_[2]);
		return false;
	}


	template <typename T> 
	bool Vector3T<T>::operator<<(const Vector& _rhs) const {
		return (data_[0] < _rhs.data_[0] && 
			data_[1] < _rhs.data_[1] && 
			data_[2] < _rhs.data_[2] ); 
	}


	template <typename T> 
	bool Vector3T<T>::operator>>(const Vector& _rhs) const {
		return (data_[0] > _rhs.data_[0] && 
			data_[1] > _rhs.data_[1] && 
			data_[2] > _rhs.data_[2] ); 
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::operator*=(const T& _s) {
			data_[0] *= _s;
			data_[1] *= _s;
			data_[2] *= _s; 
			return *this; 
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::operator/=(const T& _s) {
			assert(_s!=(T)0); 

			data_[0] /= _s;
			data_[1] /= _s;
			data_[2] /= _s; 
			return *this; 
	}


	template <typename T> 
	Vector3T<T>
		Vector3T<T>::operator*(const T& _s) const {
			return Vector(data_[0]*_s, data_[1]*_s, data_[2]*_s);
	}


	template <typename T> 
	Vector3T<T>
		Vector3T<T>::operator/(const T& _s) const {
			assert(_s!=(T)0);
			return Vector(data_[0]/_s, data_[1]/_s, data_[2]/_s);
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::operator+=(const Vector& _rhs) {
			data_[0] += _rhs.data_[0];
			data_[1] += _rhs.data_[1];
			data_[2] += _rhs.data_[2];
			return *this;
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::operator-=(const Vector& _rhs) {
			data_[0] -= _rhs.data_[0];
			data_[1] -= _rhs.data_[1];
			data_[2] -= _rhs.data_[2];
			return *this;
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::operator|=(const Vector& _rhs) {
			data_[0] *= _rhs.data_[0];
			data_[1] *= _rhs.data_[1];
			data_[2] *= _rhs.data_[2];
			return *this;
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::operator/=(const Vector& _rhs) {
			assert((_rhs.data_[0] != (T)0) &&
				(_rhs.data_[1] != (T)0) &&
				(_rhs.data_[2] != (T)0) );

			data_[0] /= _rhs.data_[0];
			data_[1] /= _rhs.data_[1];
			data_[2] /= _rhs.data_[2];

			return *this;
	}


	template <typename T> 
	Vector3T<T>
		Vector3T<T>::operator+(const Vector& _rhs) const {
			return Vector(data_[0] + _rhs.data_[0], 
				data_[1] + _rhs.data_[1], 
				data_[2] + _rhs.data_[2]);
	}


	///////////////////////////////////////////////////////////////////////////////
#define PITOANG  0.0174532925  //
	///////////////////////////////////////////////////////////////////////////////
	template <typename T>   //一个向量绕此向量旋转_angle角度
	void Vector3T<T>::rotate(Vector& _vec, const T& _angle) {
		T sinval = (T)sin(_angle*PITOANG);
		T cosval = (T)cos(_angle*PITOANG);
		T a = x * sinval;
		T b = y * sinval;
		T c = z * sinval; 
		_vec = Vector(_vec.z*b + _vec.x*cosval - _vec.y*c, 
			_vec.x*c + _vec.y*cosval - _vec.z*a,
			_vec.y*a + _vec.z*cosval - _vec.x*b);
	}
#undef PITOANG


	template <typename T> 
	void Vector3T<T>::rotatePoint(Vector& _pnt, const T& _angle) {
		T  k = x*_pnt.x + y*_pnt.y + z*_pnt.z;  //两个向量的点积

		Vector  proj (k*x, k*y, k*z);  //投影
		Vector  vect (_pnt - proj);    //旋转轴
		T       len = vect.length();   //!!

		vect.normalize();          //!!
		rotate(vect, _angle);

		vect *= len;//恢复r的长度  //!!
		_pnt  = proj + vect;
	}

	template <typename T> 
	Vector3T<T> Vector3T<T>::operator-(const Vector& _rhs) const {
		return Vector(data_[0] - _rhs.data_[0], 
			data_[1] - _rhs.data_[1], 
			data_[2] - _rhs.data_[2]);
	}


	template <typename T> 
	Vector3T<T> Vector3T<T>::operator|(const Vector& _rhs) const {
		return Vector(data_[0] * _rhs.data_[0], 
			data_[1] * _rhs.data_[1], 
			data_[2] * _rhs.data_[2]);
	}


	template <typename T> 
	Vector3T<T> Vector3T<T>::operator/(const Vector& _rhs) const {
		assert((_rhs.data_[0] != (T)0) &&
			(_rhs.data_[1] != (T)0) &&
			(_rhs.data_[2] != (T)0) );

		return Vector(data_[0] / _rhs.data_[0], 
			data_[1] / _rhs.data_[1],
			data_[2] / _rhs.data_[2]);
	}


	template <typename T> 
	Vector3T<T>
		Vector3T<T>::operator-(void) const {
			return Vector(-data_[0], -data_[1], -data_[2]);
	}


	template <typename T>                                          //dot product
	T Vector3T<T>::operator*(const Vector& _rhs) const { 
		return  data_[0] * _rhs.data_[0] +
			data_[1] * _rhs.data_[1] +
			data_[2] * _rhs.data_[2] ;
	}


	template <typename T> 	                                     //cross product
	Vector3T<T>
		Vector3T<T>::operator%(const Vector& _rhs) const {
			return Vector(data_[1]*_rhs.data_[2] - data_[2]*_rhs.data_[1],
				data_[2]*_rhs.data_[0] - data_[0]*_rhs.data_[2],
				data_[0]*_rhs.data_[1] - data_[1]*_rhs.data_[0]);
	}


	template <typename T> 
	T Vector3T<T>::norm() const {  
		return (T)
			(sqrt( data_[0]*data_[0] + data_[1]*data_[1] + data_[2]*data_[2]) );
	}


	template <typename T> 
	T Vector3T<T>::sqNorm() const {  
		return (data_[0]*data_[0] + data_[1]*data_[1] + data_[2]*data_[2]);
	}


	template <typename T>                                                  //max
	T Vector3T<T>::maxx(void) const {
		return (data_[0]>data_[1]) 
			? ((data_[0] > data_[2]) ? data_[0] : data_[2])
			: ((data_[1] > data_[2]) ? data_[1] : data_[2]);
	}


	template <typename T>                                                  //min
	T Vector3T<T>::minn(void) const {
		return (data_[0] < data_[1])
			? ((data_[0] < data_[2]) ? data_[0] : data_[2])
			: ((data_[1] < data_[2]) ? data_[1] : data_[2]);
	}


	template <typename T> 
	T Vector3T<T>::mean() const { //mean
		return  (data_[0] + data_[1] + data_[2]) / T(3); 
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::minimize(const Vector& _rhs) {
			if (_rhs.data_[0] < data_[0])  data_[0] = _rhs.data_[0];
			if (_rhs.data_[1] < data_[1])  data_[1] = _rhs.data_[1];
			if (_rhs.data_[2] < data_[2])  data_[2] = _rhs.data_[2];

			return *this;
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::maximize(const Vector& _rhs) {
			if (_rhs.data_[0] > data_[0]) data_[0] = _rhs.data_[0];
			if (_rhs.data_[1] > data_[1]) data_[1] = _rhs.data_[1];
			if (_rhs.data_[2] > data_[2]) data_[2] = _rhs.data_[2];

			return *this;
	}


	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::normalize() {  //Non int
			T len((T)(sqrt(data_[0]*data_[0] + data_[1]*data_[1] + data_[2]*data_[2])));
			if (len != (T)0) {
				len = ((T)1)/len;
				data_[0] *= len;
				data_[1] *= len;
				data_[2] *= len;
			}

			return *this;
	}


	// vectorize: set each component of vector into same value -- _s
	template <typename T> 
	Vector3T<T>&
		Vector3T<T>::vectorize(const T& _s) {
			data_[0] = data_[1] = data_[2] = _s; 
			return *this;
	} 


	///////////////////////////////////////////////////////////////////////////////
	// I/O of vector
	///////////////////////////////////////////////////////////////////////////////

	// stream output
	template <typename T>
	inline std::ostream&
		operator<<(std::ostream& _ostr, const Vector3T<T>& _v) {
			_ostr<<_v.data_[0]<<" "<<_v.data_[1]<<" "<<_v.data_[2];
			return _ostr;
	}

	// stream input, data should be comparted by 'blank'
	template <typename T>
	inline std::istream&
		operator>>(std::istream& _istr, Vector3T<T>& _v) {
			_istr>>_v.data_[0]>>_v.data_[1]>>_v.data_[2];
			return _istr;
	}


	///////////////////////////////////////////////////////////////////////////////
	// global operators for 3d vectors
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	template<typename T>                        // left scalar-product of vector
	inline Vector3T<T> 
		operator*(T _s, const Vector3T<T>& _v) {
			return (Vector3T<T>(_v) *= _s);
	}


	///////////////////////////////////////////////////////////////////////////////
	template<typename T>                          // global dot-product function
	inline T
		dot(const Vector3T<T>& _v0, const Vector3T<T>& _v1) {
			return (_v0 * _v1); 
	}


	///////////////////////////////////////////////////////////////////////////////
	template<typename T>                        // global cross-product function
	inline Vector3T<T> 
		cross(const Vector3T<T>& _v0, const Vector3T<T>& _v1) {
			return (_v0 % _v1);
	}


	///////////////////////////////////////////////////////////////////////////////
	// some specializations of 3d vector template class
	///////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////// 
	typedef Vector3T<signed char>          Vec3c;  // 3-byte signed vector
	typedef Vector3T<unsigned char>        Vec3uc; // 3-byte unsigned vector
	typedef Vector3T<signed short int>     Vec3s;  // 3-short signed vector
	typedef Vector3T<unsigned short int>   Vec3us; // 3-short unsigned vector
	typedef Vector3T<signed int>           Vec3i;  // 3-int signed vector
	typedef Vector3T<unsigned int>         Vec3ui; // 3-int unsigned vector
	typedef Vector3T<float>                Vec3f;  // 3-float vector
	typedef Vector3T<double>               Vec3d;  // 3-double vector



} //namespace 

#endif 