///////////////////////////////////////////////////////////////////////////////
#ifndef __WMQ_MESH_IO_READER_WRITER_H__
#define __WMQ_MESH_IO_READER_WRITER_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include<string>

namespace MeshN {


///////////////////////////////////////////////////////////////////////////////
template <class Mesh>
class ReaderWriterT {
public:
	typedef typename Mesh::Coord             Coord;
	typedef typename Mesh::VertexHandle      VertexHandle;
	typedef typename Mesh::FacetHandle       FacetHandle;
	typedef typename Mesh::HalfedgeHandle    HalfedgeHandle;

public:
	ReaderWriterT() : mesh_(NULL) {}
	ReaderWriterT(Mesh* _mesh) : mesh_(_mesh) {}

public:
	void set_mesh(Mesh* _ptr_mesh);    //set ptr to mesh

	bool off_reader(const char* _filename);  //reading
	//bool off_writer(char* _filename);  //write to file
	//bool ogl_writer(bool  _orient=true, bool _smooth=false);  //drawing by OGL

private:
	Mesh*  mesh_;  //mesh object
};

} //namespace

#endif