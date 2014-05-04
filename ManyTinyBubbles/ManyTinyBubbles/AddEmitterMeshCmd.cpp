#include "AddEmitterMeshCmd.h"

#include <maya/MItSelectionList.h>
#include <maya/MGlobal.h>
#include <maya/MDagPath.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>

#include "GlobalState.h"


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

AddEmitterMeshCmd::AddEmitterMeshCmd()
{
}

AddEmitterMeshCmd::~AddEmitterMeshCmd()
{
}


////////////////////////////////////////////////////
// doIt()
////////////////////////////////////////////////////
MStatus AddEmitterMeshCmd::doIt( const MArgList& args )
{
	// create a selection list to convert string to DAGPathObject
	MSelectionList selection_list;
	MGlobal::getActiveSelectionList( selection_list );

	// store the seleceted object's name in GlobalState
	MItSelectionList it( selection_list );
	while ( !it.isDone() ) {
		MDagPath dag_path ;
		MObject	component;
		it.getDagPath( dag_path, component );
		
		MFnDependencyNode fn( dag_path.node() );
		std::string obj_name = fn.name().asChar();

		// if the dag path object is a transform, get the shape
		bool found_shape = true;
		if ( dag_path.apiType() == MFn::kTransform ) {
			MStatus stat = dag_path.extendToShape();
			if ( stat != MStatus::kSuccess ) {
				found_shape = false;
			}
		}

		GlobalState::setSelectedObject( obj_name );

		// create a mesh function object if we have a shape and the node is compatible
		if ( found_shape && dag_path.hasFn( MFn::kMesh ) ) {
			// store mesh name
			GlobalState::setSelectedObject( obj_name );

			// get MFnMesh from MDagPath
			MFnMesh mesh_surface( dag_path );

			////////////////////////////////////////////////////
			// get surface data from Maya mesh object
			////////////////////////////////////////////////////

			// get number of triangles in mesh and the vertex indices for each triangle
			// triangle_vertex_indices will be three times larger than triangle_num
			MIntArray minta_tri_num;
			MIntArray minta_tri_vertex_indices;
			mesh_surface.getTriangles( minta_tri_num,
									   minta_tri_vertex_indices );
			unsigned int triangle_vertex_num = minta_tri_vertex_indices.length();

			// copy triangle_vertex_indices into a C++ vector
			std::vector<int> triangle_vertex_indices;
			triangle_vertex_indices.resize( triangle_vertex_num );
			minta_tri_vertex_indices.get( &triangle_vertex_indices[0] );


			////////////////////////////////////////////////////
			// create and store face list
			////////////////////////////////////////////////////

			std::vector<int> face_list;
			for ( unsigned int i = 0; i < triangle_vertex_num; i += 3 ) {
				face_list.push_back( triangle_vertex_indices[i] + 1 );
				face_list.push_back( triangle_vertex_indices[i + 1] + 1 );
				face_list.push_back( triangle_vertex_indices[i + 2] + 1 );
			}

			// store face list data in GlobalState
			GlobalState::storeFaceList( face_list );


			////////////////////////////////////////////////////
			// create and store vertex list
			////////////////////////////////////////////////////

			std::vector<vec3> vert_list;

			MFloatPointArray mesh_vertices;
			mesh_surface.getPoints(mesh_vertices, MSpace::kWorld);

			for ( unsigned int i = 0; i < mesh_vertices.length(); ++i ) {
				MFloatPoint vertex = mesh_vertices[i];
				vec3 new_vert( vertex[VX],
							   vertex[VY],
							   vertex[VZ] );
				vert_list.push_back( new_vert );
			}

			// store vertex list data in GlobalState
			GlobalState::storePointList( vert_list );
		}

		it.next();
	}

    return MStatus::kSuccess;
}