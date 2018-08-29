
#include "KernelT.h"

namespace MeshN {

	///////////////////////////////////////////////////////////////////////////////
	//              DEFINITIONS OF KERNEL CONVENIENT FUNCTIONS
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	// get various handles for  [\halfedge structure] 
	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename MeshN::KernelT<Items>::HalfedgeHandle 
		MeshN::KernelT<Items>::next_halfedge_handle(const HalfedgeHandle& _hh)  {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return edges_[_hh.idx()>>1].halfedges_[_hh.idx()&1].next_halfedge_handle_;
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::HalfedgeHandle
		KernelT<Items>::prev_halfedge_handle(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return edges_[_hh.idx()>>1].halfedges_[_hh.idx()&1].prev_halfedge_handle_;
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::VertexHandle
		KernelT<Items>::vertex_handle(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return edges_[_hh.idx()>>1].halfedges_[_hh.idx()&1].vertex_handle_;
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::EdgeHandle
		KernelT<Items>::edge_handle(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return EdgeHandle( _hh.idx()>>1 );
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::FacetHandle
		KernelT<Items>::facet_handle(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return edges_[_hh.idx()>>1].halfedges_[_hh.idx()&1].facet_handle_;
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::HalfedgeHandle  
		KernelT<Items>::opposite_halfedge_handle(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );

			return ( !(_hh.idx()&1) )
				? HalfedgeHandle( _hh.idx()+1 )
				: HalfedgeHandle( _hh.idx()-1 );
	}

	///////////////////////////////////////////////////////////////////////////////
	// get halfedge-handles from [\vertex \edge \facet handles] 
	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> 
	typename KernelT<Items>::HalfedgeHandle
		KernelT<Items>::halfedge_handle(const VertexHandle& _vh) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );
			return vertices_[_vh.idx()].halfedge_handle_;
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> 
	typename KernelT<Items>::HalfedgeHandle
		KernelT<Items>::halfedge_handle(const EdgeHandle& _eh, int _i) {
			assert( _eh.is_valid() );
			assert( _eh.idx() < edge_size() );
			return HalfedgeHandle ( (_eh.idx()<<1) | _i );
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> 
	typename KernelT<Items>::HalfedgeHandle
		KernelT<Items>::halfedge_handle(const FacetHandle& _fh) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );
			return facets_[_fh.idx()].halfedge_handle_;
	}

	///////////////////////////////////////////////////////////////////////////////
	// rotating halfedge-handles 
	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> 
	typename KernelT<Items>::HalfedgeHandle
		KernelT<Items>::ccw_rotated(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return prev_halfedge_handle( opposite_halfedge_handle(_hh) );
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> 
	typename KernelT<Items>::HalfedgeHandle
		KernelT<Items>::cw_rotated(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return opposite_halfedge_handle( next_halfedge_handle(_hh) );
	}

	///////////////////////////////////////////////////////////////////////////////
	// handles <--> refereces
	/////////////////////////////////////////////////////////////////////////////// 
	template  <class Items>
	typename KernelT<Items>::Halfedge&
		KernelT<Items>::halfedge_ref(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );

			int idx (_hh.idx()); 
			return edges_[_hh.idx()>>1].halfedges_[_hh.idx()&1];
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::Vertex&
		KernelT<Items>::vertex_ref(const VertexHandle& _vh) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );
			return vertices_[_vh.idx()];
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::Edge&
		KernelT<Items>::edge_ref(const EdgeHandle& _eh) {
			assert( _eh.is_valid() );
			assert( _eh.idx() < edge_size() );
			return edges_[_eh.idx()];
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	typename KernelT<Items>::Facet&
		KernelT<Items>::facet_ref(const FacetHandle& _fh) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );
			return facets_[_fh.idx()];
	}


	///////////////////////////////////////////////////////////////////////////////
	//               DEFINITIONS OF COMMON MEMBER FUNCTIONS
	///////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////// 
	// clear and allign
	/////////////////////////////////////////////////////////////////////////////// 
	template <class Items>                                                 // clear
	void KernelT<Items>::clear(void) {
		vertices_.clear();
		edges_.clear();
		facets_.clear(); 
	}

	/////////////////////////////////////////////////////////////////////////////// 
	template <class Items>                                                // assign
	typename KernelT<Items>::This&  
		KernelT<Items>::operator=(const This& _rhs) { 
			vertices_ = _rhs.vertices_; 
			edges_    = _rhs.edges_; 
			facets_   = _rhs.facets_; 

			return *this; 
	}

	///////////////////////////////////////////////////////////////////////////////
	// is boundary ?
	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> bool
		KernelT<Items>::is_boundary_halfedge(const HalfedgeHandle& _hh )  {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return !facet_handle(_hh).is_valid();
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> bool
		KernelT<Items>::is_boundary_edge(const EdgeHandle& _eh)  {
			assert( _eh.is_valid() );
			assert( _eh.idx() < edge_size() );
			return (is_boundary_halfedge( halfedge_handle(_eh, 0)) ||
				is_boundary_halfedge( halfedge_handle(_eh, 1)) );
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> bool
		KernelT<Items>::is_boundary_vertex(const VertexHandle& _vh) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );

			HalfedgeHandle hh (halfedge_handle(_vh)); 

			if ( !hh.is_valid() ) return true; // found: isolated point 

			HalfedgeHandle cursor (hh);
			do {
				if ( is_boundary_halfedge(cursor) )  return true;
				cursor = cw_rotated( cursor );
			} while ( cursor != hh ); 

			return false;
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> bool
		KernelT<Items>::is_boundary_facet(const FacetHandle& _fh) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );

			HalfedgeHandle hh (halfedge_handle(_fh));
			assert( hh.is_valid() );

			HalfedgeHandle cursor (hh);
			do {
				if ( is_boundary_halfedge( opposite_halfedge_handle(cursor)) )
					return true;
				cursor = next_halfedge_handle( cursor );
			} while ( cursor != hh );

			return false;
	}

	///////////////////////////////////////////////////////////////////////////////
	// numbers of vertices/facets/edges in mesh
	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	int KernelT<Items>::vertex_number(void) {
		int            num(0);
		VertexIterator vi(vertices_.begin());

		for ( ; vi!=vertices_.end(); ++vi ) {
			if ( vi->status_.is_deleted() ) continue;
			if ( vi->status_.is_hidden() )  continue;
			++num;
		}

		return num;
	}

	////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	int KernelT<Items>::edge_number(void) {
		int          num(0);
		EdgeIterator ei(edges_.begin());

		for ( ; ei!=edges_.end(); ++ei ) {
			if ( ei->status_.is_deleted() ) continue;
			if ( ei->status_.is_hidden() )  continue;
			++num;
		}

		return num;
	}

	////////////////////////////////////////////////////////////////////////////
	template  <class Items>
	int KernelT<Items>::facet_number(void) {
		int num(0);
		FacetIterator fi(facets_.begin());

		for (; fi!=facets_.end(); ++fi ) {
			if ( fi->status_.is_deleted() ) continue;
			if ( fi->status_.is_hidden() )  continue;
			++num;
		}

		return num;
	}

	///////////////////////////////////////////////////////////////////////////////
	template <class Items>                            // get coord from vertex 
	typename KernelT<Items>::Coord& 
		KernelT<Items>::coord(const VertexHandle& _vh) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );
			return vertex_ref(_vh).coord_;
	} 



	///////////////////////////////////////////////////////////////////////////////
	//           DEFINITIONS OF MANAGER FUNCIONS OF STATUS BITS
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	// bit - deleted
	/////////////////////////////////////////////////////////////////////////////// 
	template  <class Items> void
		KernelT<Items>::set_deleted(const HalfedgeHandle& _hh, bool _b) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			halfedge_ref(_hh).status_.set_deleted(_b);
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> void
		KernelT<Items>::set_deleted(const VertexHandle& _vh, bool _b) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );
			vertex_ref(_vh).status_.set_deleted(_b);
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> void
		KernelT<Items>::set_deleted(const EdgeHandle& _eh, bool _b) {
			assert( _eh.is_valid() );
			assert( _eh.idx() < edge_size() );

			edge_ref(_eh).status_.set_deleted(_b);
			edge_ref(_eh).halfedges_[0].status_.set_deleted(_b);
			edge_ref(_eh).halfedges_[1].status_.set_deleted(_b);
	}

	///////////////////////////////////////////////////////////////////////////////
	template <class Items> void
		KernelT<Items>::set_deleted(const FacetHandle& _fh, bool _b) {
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );
			facet_ref(_fh).status_.set_deleted(_b);
	}

	///////////////////////////////////////////////////////////////////////////////
	template <class Items> bool
		KernelT<Items>::is_deleted(const VertexHandle& _vh) {
			assert( _vh.is_valid() );
			assert( _vh.idx() < vertex_size() );
			return vertex_ref(_vh).status_.is_deleted();
	}

	///////////////////////////////////////////////////////////////////////////////
	template <class Items> bool
		KernelT<Items>::is_deleted(const HalfedgeHandle& _hh) {
			assert( _hh.is_valid() );
			assert( _hh.idx() < halfedge_size() );
			return halfedge_ref(_hh).status_.is_deleted();
	}

	///////////////////////////////////////////////////////////////////////////////
	template <class Items> bool
		KernelT<Items>::is_deleted(const EdgeHandle& _eh) {
			assert( _eh.is_valid() );
			assert( _eh.idx() < edge_size() );
			return edge_ref(_eh).status_.is_deleted();
	}

	///////////////////////////////////////////////////////////////////////////////
	template <class Items> bool
		KernelT<Items>::is_deleted(const FacetHandle& _fh) { 
			assert( _fh.is_valid() );
			assert( _fh.idx() < facet_size() );
			return facet_ref(_fh).status_.is_deleted();
	}

	///////////////////////////////////////////////////////////////////////////////
	// bit - hidden
	///////////////////////////////////////////////////////////////////////////////


	//............


	///////////////////////////////////////////////////////////////////////////////
	//              DEFINITIONS OF BUILDING FUNCTIONS OF MESH
	//        TO CONSTRUCT HALFEDGE SYSTEM IN MEMORY FROM FILE DATA
	/////////////////////////////////////////////////////////////////////////////// 
	template <class Items>                                            // add Vertex
	typename KernelT<Items>::VertexHandle
		KernelT<Items>::add_vertex(void) {
			vertices_.push_back( Vertex() ); 
			return VertexHandle( vertices_.size()-1 );
	} 

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>                                           // add Vertex
	typename KernelT<Items>::VertexHandle
		KernelT<Items>::add_vertex(const Coord& _coord) {
			Vertex vertex; 
			vertex.coord_ = _coord; 
			vertex.oriCoord_ = _coord; //lzhu
			vertex.sumNormal_ = MathN::Vec3f(0, 0, 0);
			vertex.weights_ = MathN::Vec3f(0, 0, 0);
			vertices_.push_back( vertex ); 

			return VertexHandle( vertices_.size()-1 );
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>                                    // add facet to mesh
	typename KernelT<Items>::FacetHandle
		KernelT<Items>::add_facet(const std::vector<VertexHandle>& _vhs) { 
			std::vector<HalfedgeHandle>   halfedge_handles; //记录新面片的半边号!
			std::vector<bool>             flags;            //上述半边已在网格中?
			FacetHandle                   new_facet_handle (facets_.size());  //新面柄 
			int  i, j, cnt(_vhs.size());                    //循环变量&循环记数器

			assert( cnt >= 3 );  //无顶点数少于3的面!

			///验证所有顶点是否都为边界顶点///////
			for ( i=0; i<cnt; ++i ) {
				//assert( _vhs[i].is_valid() );
				if ( !is_boundary_vertex(_vhs[i]) ) {
					std::cout<<"Error: Complex vertex ["<<_vhs[i]<<"].\n";
					return FacetHandle();
					////exit(0);
				}
			} 

			///注册所有的半边,并加入缺失的半边(与边)
			for (j=cnt-1, i=0; i<cnt; j=(j+1)%cnt, ++i) { 
				HalfedgeHandle heh = _find_halfedge_handle(_vhs[j], _vhs[i]); 

				if ( heh.is_valid() ) {  
					//如果半边已在网格中,则验证其是否为边界半边
					if ( ! is_boundary_halfedge(heh) ) { //为复杂边
						std::cout<<"Error: Complex edge: <";
						std::cout<<_vhs[j]<<", "<<_vhs[i]<<"("<<heh<<")"<<"> !\n";
						exit(0);
					}
					halfedge_handles.push_back( heh );  //注册之
					flags.push_back(true); //标记之
				}  
				else {
					//半边不在网格中,需要加入缺失的边(半边), 并注册之 
					halfedge_handles.push_back(
						halfedge_handle( add_edge(_vhs[j], _vhs[i]), 1) );
					//注意:上面加入边的缺失半边为第二条半边!!!!!!!!!!
					flags.push_back(false);
				}

				//置所有半边的面片号
				halfedge_ref(halfedge_handles[i]).facet_handle_ = new_facet_handle;
			}

			///建立与新添加面片相关的半边拓扑关系
			for (i=0, j=1; i<cnt; ++i, j=(j+1)%cnt) {
				int flag(0);  //标记相继两条半边的存在(新/旧)关系

				if ( !flags[i] )  flag |= 1; //相继半边的关系
				if ( !flags[j] )  flag |= 2;

				Halfedge& he0 = halfedge_ref( halfedge_handles[i] );
				switch (flag){
				case 0: { //两条半边均不是刚刚生成的
					//如果位置错乱, 则调整之
					if ( next_halfedge_handle(halfedge_handles[i]) != halfedge_handles[j] )
						_fix_halfedge_order( halfedge_handles[i], halfedge_handles[j] );
					break;
						} //case 0

				case 1: { //前一条半边新的,后一条是老的
					he0.next_halfedge_handle_ = halfedge_handles[j];

					halfedge_ref(opposite_halfedge_handle(halfedge_handles[i])).
						prev_halfedge_handle_ = prev_halfedge_handle( halfedge_handles[j] );

					halfedge_ref(prev_halfedge_handle(halfedge_handles[j])).
						next_halfedge_handle_ = opposite_halfedge_handle( halfedge_handles[i] );

					halfedge_ref(halfedge_handles[j]).
						prev_halfedge_handle_ = halfedge_handles[i];

					//NOTE:顶点(_vns[i])的第一条(边界)半边保持不变!!
					break;
						} //case 1

				case 2: {//前一条半边是老的, 而后一条则是新近生成的
					halfedge_ref(halfedge_handles[j]).
						prev_halfedge_handle_ = halfedge_handles[i];

					halfedge_ref(opposite_halfedge_handle(halfedge_handles[j])).
						next_halfedge_handle_ = next_halfedge_handle( halfedge_handles[i] );

					halfedge_ref(next_halfedge_handle(halfedge_handles[i])).
						prev_halfedge_handle_ = opposite_halfedge_handle( halfedge_handles[j] );

					he0.next_halfedge_handle_ = halfedge_handles[j];

					//NOTE: 半边_opposite_halfedge_handle(halfedge_handles[j])是边界半边!! 
					vertex_ref(_vhs[i]).halfedge_handle_ =
						opposite_halfedge_handle(halfedge_handles[j]);

					break;
						} //case 2

				case 3: {//两条半边均是新近产生的
					he0.next_halfedge_handle_ = halfedge_handles[j];
					halfedge_ref(halfedge_handles[j]).
						prev_halfedge_handle_ = halfedge_handles[i];

					if ( !halfedge_handle(_vhs[i]).is_valid() ) { //???????
						//如果(_vns[i])是孤立顶点
						halfedge_ref(opposite_halfedge_handle(halfedge_handles[i])).
							prev_halfedge_handle_ = opposite_halfedge_handle(halfedge_handles[j]);

						halfedge_ref(opposite_halfedge_handle(halfedge_handles[j])).
							next_halfedge_handle_ = opposite_halfedge_handle(halfedge_handles[i]);

						//Note: 半边_opposite_halfedge_handle(halfedge_handles[j])为边界半边
						vertex_ref(_vhs[i]).halfedge_handle_
							= opposite_halfedge_handle( halfedge_handles[j] );  

					} else {
						//否则,顶点(_vns[i])不是孤立顶点,它存在着其它连接关系
						HalfedgeHandle& v_heh = halfedge_handle( _vhs[i] ); 

						//Note: 这时半边[_halfedgeHandle(_vhs[i])]应为边界半边
						assert( is_boundary_halfedge(v_heh) );

						halfedge_ref(opposite_halfedge_handle(halfedge_handles[j])).
							next_halfedge_handle_ = next_halfedge_handle( v_heh );

						halfedge_ref(next_halfedge_handle(v_heh)).
							prev_halfedge_handle_ = opposite_halfedge_handle(halfedge_handles[j]);

						halfedge_ref(v_heh).next_halfedge_handle_ =
							opposite_halfedge_handle( halfedge_handles[i] );

						halfedge_ref(opposite_halfedge_handle(halfedge_handles[i])).
							prev_halfedge_handle_ = v_heh;

						//NOTE:这时顶点的第一半边仍为边界半边,无需调整
					}//case 3
						}//switch
				}//end of for

				_adjust_to_boundary_halfedge( _vhs[i] );
			}

			//产生新面片,注意其中的半边拓扑关系没有完全确定下来
			Facet facet;
			facet.halfedge_handle_ = halfedge_handles[0];
			facets_.push_back( facet ); 

			return new_facet_handle;
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>                              // add a tri-facet to mesh
	typename KernelT<Items>::FacetHandle
		KernelT<Items>::add_facet(const VertexHandle& _vh0, ////&
		const VertexHandle& _vh1,
		const VertexHandle& _vh2 )
	{ 
		VertexHandle vhs[3] = { _vh0, _vh1, _vh2 };

		///验证所有顶点是否都为边界顶点//////////////////////
		//assert( vhs[0].is_valid() && vhs[1].is_valid() && vhs[2].is_valid() );
		assert( is_boundary_vertex(vhs[0]) && is_boundary_vertex(vhs[1]) &&
			is_boundary_vertex(vhs[2]) );

		HalfedgeHandle  heHandles[3]; //记录新面片的半边号!
		bool            flags[3];     //上述半边已在网格中?
		FacetHandle     newFacetHandle ( facets_.size() );  //被添入面片的序号!
		int             i, j;     //循环变量&循环记数器

		///注册所有的半边,并加入缺失的半边(与边)/////////////
		for (j=2, i=0; i<3; j=(j+1)%3, ++i) {
			HalfedgeHandle heh = _find_halfedge_handle(vhs[j], vhs[i]); 
			//如果半边已在网格中,验证其是否为边界半边
			if ( heh.is_valid() ) { 
				assert( is_boundary_halfedge(heh) ); //不为复杂边
				heHandles[i] = heh;  //注册之
				flags[i]     = true; //标记之
			}
			//半边不在网格中,需要加入缺失的边(半边), 并注册之
			else {
				heHandles[i] = halfedge_handle( add_edge(vhs[j], vhs[i]), 1);
				//注意:上面加入边的第一条半边为缺失半边,而非第二条半边!!!!
				flags[i] = false;
			}
			//置所有半边的面片号
			halfedge_ref(heHandles[i]).facet_handle_ = newFacetHandle;
		}

		///建立与新添加面片相关的半边拓扑关系////////////////
		for (i=0, j=1; i<3; ++i, j=(j+1)%3) {
			int   flag(0);  //标记相继两条半边的存在(新/旧)关系
			Halfedge& he0 = halfedge_ref( heHandles[i] );

			if (!flags[i])  flag |= 1; //相继半边的关系
			if (!flags[j])  flag |= 2;

			switch (flag) 	{
			case 0: { //两条半边均不是刚刚生成的
				//如果位置错乱,调整之
				if ( next_halfedge_handle(heHandles[i]) != heHandles[j] ) {
					_fix_halfedge_order(heHandles[i], heHandles[j]);
				}
				break;
					}//case 0

			case 1: {//前一条半边是新的,后一条是老的
				he0.next_halfedge_handle_ = heHandles[j];

				halfedge_ref(opposite_halfedge_handle(heHandles[i])).prev_halfedge_handle_
					= prev_halfedge_handle(heHandles[j]);

				halfedge_ref(prev_halfedge_handle(heHandles[j])).next_halfedge_handle_
					= opposite_halfedge_handle(heHandles[i]);

				halfedge_ref(heHandles[j]).prev_halfedge_handle_ = heHandles[i];

				//NOTE:顶点(vhs[i])的第一条(边界)半边保持不变!!
				break;
					}//case 1

			case 2: {//前一条半边是老的, 而后一条则是新近生成的
				halfedge_ref(heHandles[j]).prev_halfedge_handle_ = heHandles[i];

				halfedge_ref(opposite_halfedge_handle(heHandles[j])).
					next_halfedge_handle_ = next_halfedge_handle( heHandles[i] );

				halfedge_ref(next_halfedge_handle(heHandles[i])).
					prev_halfedge_handle_ = opposite_halfedge_handle( heHandles[j] );

				he0.next_halfedge_handle_ = heHandles[j];

				//NOTE: 半边[_opposite_halfedge_handle(heHandles[j])]是边界半边!!
				vertex_ref(vhs[i]).halfedge_handle_ 
					= opposite_halfedge_handle( heHandles[j] );
				break;
					}//case 2

			case 3: {//两条半边均是新近产生的
				he0.next_halfedge_handle_ = heHandles[j];
				halfedge_ref(heHandles[j]).prev_halfedge_handle_ = heHandles[i];

				if ( !halfedge_handle(vhs[i]).is_valid() ) { 
					//如果(_vns[i])是孤立顶点
					halfedge_ref(opposite_halfedge_handle(heHandles[i])).
						prev_halfedge_handle_ = opposite_halfedge_handle( heHandles[j] );

					halfedge_ref(opposite_halfedge_handle(heHandles[j])).
						next_halfedge_handle_ = opposite_halfedge_handle( heHandles[i] );

					//Note: 半边[_opposite_halfedge_handle(heHandles[j])]为边界半边
					vertex_ref(vhs[i]).halfedge_handle_ 
						= opposite_halfedge_handle( heHandles[j] ); 
				} else {
					//否则,顶点(_vns[i])不是孤立顶点,它存在着其它连接关系
					HalfedgeHandle& v_hen = halfedge_handle( vhs[i] ); 

					//Note: 这时半边[vertex_HalfedgeHandle(_vns[i])]应为边界半边
					assert( is_boundary_halfedge(v_hen) );

					halfedge_ref(opposite_halfedge_handle(heHandles[j])).
						next_halfedge_handle_ = next_halfedge_handle( v_hen );

					halfedge_ref(next_halfedge_handle(v_hen)).
						prev_halfedge_handle_ = opposite_halfedge_handle(heHandles[j]);

					halfedge_ref(v_hen).next_halfedge_handle_
						= opposite_halfedge_handle(heHandles[i]);

					halfedge_ref(opposite_halfedge_handle(heHandles[i])).
						prev_halfedge_handle_ = v_hen;

					//NOTE:这时顶点的第一半边仍为边界半边,无需调整
				}//case 3
					}//switch
			}//end of for

			_adjust_to_boundary_halfedge(vhs[i]);
		} 

		//产生新面片,注意其中的半边拓扑关系没有完全确定下来
		Facet facet;
		facet.halfedge_handle_ = heHandles[0];
		facets_.push_back( facet );  

		return newFacetHandle; 
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> 
	typename KernelT<Items>::EdgeHandle
		KernelT<Items>::add_edge(const VertexHandle& _vh0, const VertexHandle& _vh1) {  
			Edge edge; 

			edge.halfedges_[0].vertex_handle_ = _vh0; 
			edge.halfedges_[1].vertex_handle_ = _vh1; 
			edges_.push_back(edge); 

			return EdgeHandle( edges_.size()-1 );
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items> 
	typename KernelT<Items>::HalfedgeHandle
		KernelT<Items>::
		_find_halfedge_handle(const VertexHandle& _begin, const VertexHandle& _end) {
			assert( _begin.idx()<vertices_.size() && _end.idx()<vertices_.size() );

			HalfedgeHandle hh (halfedge_handle(_end));
			if ( !hh.is_valid() )  //heh1 is a isolated point handle
				return HalfedgeHandle();

			HalfedgeHandle cursor (hh); 		
			do {
				if ( vertex_handle(opposite_halfedge_handle(cursor)) == _begin )
					return cursor;
				cursor = opposite_halfedge_handle( next_halfedge_handle(cursor) );
			} while ( cursor != hh );

			return HalfedgeHandle();
	}

	///////////////////////////////////////////////////////////////////////////////
	template  <class Items>  void
		typename KernelT<Items>::
		_fix_halfedge_order(const HalfedgeHandle& _hh0, const HalfedgeHandle& _hh1) {
			//将被错误放置的两条相继半边调整为正常邻接关系;  参数: 两条相继的半边号
			HalfedgeHandle  fore_begin,          fore_end(_hh0);
			HalfedgeHandle  this_begin(_hh1),    this_end;
			HalfedgeHandle  back_begin(next_halfedge_handle(_hh0));
			HalfedgeHandle  back_end  (prev_halfedge_handle(_hh1));

			//求this_end和fore_begin
			this_end = opposite_halfedge_handle( this_begin );
			while ( !is_boundary_halfedge( this_end ) ) {
				this_end = opposite_halfedge_handle( next_halfedge_handle(this_end) );
				if ( this_end == fore_end ) {
					std::cout<<"Error: Vertex's topology wrong !";
					return; //exit(0);
				}
			};

			fore_begin = next_halfedge_handle( this_end );

			//调整半边关系: 首先摘除中间面片
			halfedge_ref(back_end).next_halfedge_handle_   = fore_begin;
			halfedge_ref(fore_begin).prev_halfedge_handle_ = back_end;

			//再将摘除的中间面片插入到正确位置
			halfedge_ref(fore_end).next_halfedge_handle_   = this_begin;
			halfedge_ref(this_begin).prev_halfedge_handle_ = fore_end;
			halfedge_ref(this_end).next_halfedge_handle_   = back_begin;
			halfedge_ref(back_begin).prev_halfedge_handle_ = this_end;

			vertex_ref(vertex_handle(_hh0)).halfedge_handle_ = back_end;
	}

	////////////////////////////////////////////////////////////////////////////
	//尽可能将顶点的第一条半边调整为边界入半边
	template  <class Items> 
	bool KernelT<Items>::
		_adjust_to_boundary_halfedge(const VertexHandle& _vh) {
			HalfedgeHandle hh (halfedge_handle(_vh));

			if ( is_boundary_halfedge(hh) ) return true;

			HalfedgeHandle cursor (opposite_halfedge_handle(
				next_halfedge_handle(hh)));
			do {
				if ( is_boundary_halfedge(cursor) ) {
					vertex_ref(_vh).halfedge_handle_ = cursor;
					return true;
				}
				cursor = opposite_halfedge_handle(next_halfedge_handle(cursor));
			} while ( cursor != hh );

			return false;
	}

	template<class Items>
	bool KernelT<Items>::calcMeshBox(){

		Coord min, max;

		VertexIterator vit = vertex_begin();
		min = max = vit->coord_;
		for(++vit; vit != vertex_end(); vit++){

			Scalar x = vit->coord_[0], y = vit->coord_[1], z = vit->coord_[2];

			if(x < min[0]) min[0] = x;
			else if(x > max[0]) max[0] = x;

			if(y < min[1]) min[1] = y;
			else if(y > max[2]) max[2] = y;

			if(z < min[2]) min[2] = z;
			else if(z > max[2]) max[2] = z;


		}

		setMeshBox(min,max);

		return true;
	}


} //namespace MeshN


