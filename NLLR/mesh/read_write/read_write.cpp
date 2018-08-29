///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include<string>
#include "read_write.h"
namespace MeshN {
	///////////////////////////////////////////////////////////////////////////////
	// implementation
	///////////////////////////////////////////////////////////////////////////////
	template <class Mesh>
	void ReaderWriterT<Mesh>::set_mesh(Mesh* _ptr_mesh) {
		mesh_ = _ptr_mesh;
	}

	///////////////////////////////////////////////////////////////////////////////
	template <class Mesh>
	bool ReaderWriterT<Mesh>::off_reader(const char* _filename) { 
		if ( mesh_ == NULL ) { 
			std::cout<<"\nPointer to mesh is NOT defined!\n";
			exit(0);
		}

		FILE *modelfile;

		if ((modelfile = fopen( _filename, "r") ) == NULL) {
			std::cout<<"\nError: Can't open the model file!\n";
			return false;  
		}

		char letter1 = ' ',letter2=' ',letter3=' ';
		for ( ; ; ) {             // skip over comments and blank 
			if (feof(modelfile)) { // no data in file
				std::cout<<"\nError: No data in model file!\n";
				return false;  
			}

			fscanf(modelfile, "%c %c %c", &letter1,&letter2,&letter3);
			if ((letter1 == 'O')&&(letter2=='F')&&(letter3=='F')) break;  //OK!!
		} 

		int no_vertex,no_facet,no_edge;
		//std::cout<<no_vertex;

		fscanf(modelfile,"%d %d %d",&no_vertex,&no_facet,&no_edge);
		//std::cout<<no_vertex;


		float  x,  y,  z;   // components of coordinats 
		int    p0, p1, p2,p3;  // index of points
		for(int i=0;i<no_vertex;i++)
		{
			fscanf(modelfile, "%f %f %f", &x, &y, &z);
			mesh_->add_vertex( Coord(x, y, z) );
			
		}

		for(int i=0;i<no_facet;i++)
		{ 
			fscanf(modelfile, "%d %d %d %d", &p0, &p1, &p2,&p3); //// ????????
			//mesh_->add_facet( VertexHandle(p0-1), VertexHandle(p1-1), VertexHandle(p2-1));
			mesh_->add_facet( VertexHandle(p1), VertexHandle(p2), VertexHandle(p3));
		}  // end of while

		fclose(modelfile);  

		return true;
	}


	/////////////////////////////////////////////////////////////////////////////// 
	//template <class Mesh>
	//bool SmfReaderWriterT<Mesh>::ogl_writer(bool _orient, bool _smooth){
	//	HalfedgeHandle       cshh;
	//	Mesh::FacetIterator  fit (mesh_->facet_begin()); 

	//	if (_smooth) {
	//		glShadeModel(GL_FLAT);
	//		//glDisable(GL_LIGHTING);
	//		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//		int orient =true;// (_orient) ? 1 : -1; 
	//		//mesh_->update_normals();  //mesh_->update_vertex_normals();//

	//		for ( ; fit != mesh_->facet_end(); ++fit ) {	

	//			if ( (*fit).status_.is_deleted() ) continue;  

	//			cshh = fit->halfedge_handle_;
	//			FacetHandle fh = mesh_->facet_handle(cshh);
	//			//const VertexHandle& vh0 = mesh_->vertex_handle(cshh);
	//			glNormal3fv( mesh_->normal(fh)*orient );
	//			glBegin(GL_TRIANGLES); 
	//			do {
	//				const VertexHandle& vh = mesh_->vertex_handle(cshh);

	//				//glNormal3fv( mesh_->normal(vh)*orient );
	//				//glColor3fv(mesh_->normal(vh));
	//				glVertex3fv( mesh_->coord(vh) );
	//				cshh = mesh_->next_halfedge_handle(cshh);
	//			} while ( cshh != fit->halfedge_handle_ );
	//			glEnd();//*/
	//		}

	//		glDisable(GL_LIGHTING);
	//		for(fit = mesh_->facet_begin(); fit != mesh_->facet_end(); fit++){

	//			cshh = fit->halfedge_handle_;
	//			glLineWidth(3);
	//			glBegin(GL_LINES); //??????
	//			do {
	//				const VertexHandle& vh = mesh_->vertex_handle(cshh);

	//				glColor3f(0.0,0,1.0);
	//				glVertex3fv( mesh_->coord(vh) );
	//				glVertex3fv( mesh_->normal(vh)*0.15 + mesh_->coord(vh) );
	//				cshh = mesh_->next_halfedge_handle(cshh);
	//			} while ( cshh != fit->halfedge_handle_ );
	//			glEnd();//
	//		}
	//		glEnable(GL_LIGHTING);

	//	} else {
	//		for ( ; fit!=mesh_->facet_end(); ++fit ) {	
	//			if ( (*fit).status_.is_deleted() ) continue; 

	//			// calculating normals for facets
	//			const HalfedgeHandle& hh0 = fit->halfedge_handle_;
	//			const HalfedgeHandle& hh1 = mesh_->next_halfedge_handle(hh0);
	//			const HalfedgeHandle& hh2 = mesh_->next_halfedge_handle(hh1);

	//			const Coord& cd0 = mesh_->coord( mesh_->vertex_handle(hh0) );
	//			const Coord& cd1 = mesh_->coord( mesh_->vertex_handle(hh1) );
	//			const Coord& cd2 = mesh_->coord( mesh_->vertex_handle(hh2) );  

	//			Coord nm0 = (cd0-cd1) % (cd0-cd2);

	//			cshh = hh0;
	//			glBegin(GL_POLYGON); 
	//			(_orient) ? glNormal3fv(-nm0.normalize())
	//				: glNormal3fv(nm0.normalize());
	//			do {
	//				glVertex3fv( mesh_->coord(mesh_->vertex_handle(cshh)) );
	//				cshh = mesh_->next_halfedge_handle(cshh);
	//			} while ( cshh != hh0 );
	//			glEnd();
	//		}
	//	} 

	//	return true;
	//}

	
} //namespace
