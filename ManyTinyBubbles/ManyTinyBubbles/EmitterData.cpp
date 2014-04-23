#include "EmitterData.h"

#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>

#include "SDFGen/makelevelset3.h"
#include "Convenience.h"


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

EmitterData::EmitterData()
{
}

EmitterData::~EmitterData()
{
}


////////////////////////////////////////////////////
// set m_name
////////////////////////////////////////////////////
void EmitterData::init( MString name )
{
	m_name = name;
}


////////////////////////////////////////////////////
// generate list of positions on mesh where bubbles can emit from
// basically, fill m_source_pos_list
////////////////////////////////////////////////////
void EmitterData::createEmissionPositionsOnMesh( const unsigned int& voxel_num,
												 const unsigned int& melting_rate )
{
	// TODO: see if we can remove some of the data copying we're currently doing
	// TODO: fill m_source_pos_list directly instead of filling source_pos_list and copying
	// TODO: break this method up into multiple smaller methods b/c it's too long and difficult to work in

	// list of possible bubble generation locations on mesh
	std::vector<vec3> source_pos_list;


	////////////////////////////////////////////////////
	// first, get Maya mesh object from mesh name
	////////////////////////////////////////////////////

	// select Maya object by name
	MSelectionList maya_sel_list;
	MGlobal::getSelectionListByName( m_name, maya_sel_list );

	// get path to a mesh DAG node
	MDagPath dag_node_path;
	maya_sel_list.getDagPath( 0, dag_node_path );

	// if dag_node_path points to a transform node instead of a shape node, try to extend the DAG path to reach a shape node
	bool node_has_shape = true;
	if ( dag_node_path.apiType() == MFn::kTransform ) {
		MStatus stat = dag_node_path.extendToShape();
		if ( stat != MStatus::kSuccess ) {
			node_has_shape = false;
		}
	}

	// if we were able to find the shape node for the input mesh, and the node supports the kMesh function set
	if ( node_has_shape && dag_node_path.hasFn( MFn::kMesh ) ) {

		// get MFnMesh from MDagPath
		MFnMesh mesh_surface( dag_node_path );


		////////////////////////////////////////////////////
		// next, get surface data from Maya mesh object
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

		// create min and max points for bounding box
		// min_bounds is the minimum x, y, z of all vertices
		// max_bounds is the maximum x, y, z of all vertices
		// space enclosed between min_bounds and max_bounds enclose all vertices
		Vec3f min_bounds( std::numeric_limits<float>::max(),
						  std::numeric_limits<float>::max(),
						  std::numeric_limits<float>::max() );
		Vec3f max_bounds( -std::numeric_limits<float>::max(),
						  -std::numeric_limits<float>::max(),
						  -std::numeric_limits<float>::max() );

		
		////////////////////////////////////////////////////
		// create face list
		// use Vec3ui for compatibility with SDFGen level set library
		////////////////////////////////////////////////////

		std::vector<Vec3ui> face_list;

		// create face list
		for ( unsigned int i = 0; i < triangle_vertex_num; i += 3 ) {
			face_list.push_back( Vec3ui( triangle_vertex_indices[i],
										 triangle_vertex_indices[i + 1],
										 triangle_vertex_indices[i + 2] ) );
		}


		////////////////////////////////////////////////////
		// create vertex list
		// use Vec3f for compatibility with SDFGen level set library
		////////////////////////////////////////////////////

		std::vector<Vec3f> vert_list;

		MFloatPointArray mesh_vertices;
		mesh_surface.getPoints(mesh_vertices, MSpace::kWorld);

		for ( unsigned int i = 0; i < mesh_vertices.length(); ++i ) {
			MFloatPoint vertex = mesh_vertices[i];
			Vec3f new_vert( vertex[VX],
							vertex[VY],
							vertex[VZ] );
			vert_list.push_back( new_vert );

			// update_minmax is a method in vecL.h of SDFGen library
			// update bounding box points
			update_minmax( new_vert, min_bounds, max_bounds );
		}


		////////////////////////////////////////////////////
		// bounding box stuff
		////////////////////////////////////////////////////

		float bounding_box_size_x = max_bounds[VX] - min_bounds[VX];
		float bounding_box_size_y = max_bounds[VY] - min_bounds[VY];
		float bounding_box_size_z = max_bounds[VZ] - min_bounds[VZ];

		// get longest bounding box side
		float bounding_box_max_length = bounding_box_size_x;
		if ( bounding_box_size_y > bounding_box_max_length ) {
			bounding_box_max_length = bounding_box_size_y;
		}
		if ( bounding_box_size_z > bounding_box_max_length ) {
			bounding_box_max_length = bounding_box_size_z;
		}

		// TODO: ask about dx computation

		// compute reasonable dx for computing level sets
		// dx represents the level set resolution
		float bounding_box_volume = bounding_box_size_x * bounding_box_size_y * bounding_box_size_z;
		float voxel_volume = bounding_box_volume / voxel_num; // voxel_num is a method argument
		float dx = pow( voxel_volume, ( 1.0f / 3.0f ) );

		// TODO: ask why we're padding the bounding box

		// add padding around bounding box
		Vec3f unit( 1.0f, 1.0f, 1.0f );
		min_bounds -= 0.001f * unit;
		max_bounds += 0.001f * unit;


		////////////////////////////////////////////////////
		// compute signed distance function
		////////////////////////////////////////////////////

		// TODO: WTF?!
		// number of sample points per x, y, z direction?
		Vec3ui sizes = Vec3ui( ( max_bounds - min_bounds ) / dx ) + Vec3ui( 2, 2, 2 );

		// Array3f is a data structure defined in array3.h of SDFGen library
		Array3f signed_distance_values;

		// make_level_set3 is a method in makelevelset3.cpp of SDFGen library
		make_level_set3( face_list,
						 vert_list,
						 min_bounds,
						 dx,
						 sizes[VX],
						 sizes[VY],
						 sizes[VZ],
						 signed_distance_values );

		// create bubble generation positions from signed distance function comuted in make_level_set3
		// signed_distance_values.a is an ArrayT
		for ( unsigned int i = 0; i < signed_distance_values.a.size(); ++i ) {

			// bubble generation position will be on outer surface of mesh or within some threshold above it
			if ( signed_distance_values.a[i] >= 0.0f && signed_distance_values.a[i] < dx / 2.0f ) {

				// TODO: ask what ni, nj, nk are for
				// TODO: ask about this position computation logic

				int z = i / signed_distance_values.ni / signed_distance_values.nj;
				int y = ( i - z * signed_distance_values.ni * signed_distance_values.nj ) / signed_distance_values.ni;
				int x = i - z * signed_distance_values.ni * signed_distance_values.nj - y * signed_distance_values.ni;

				// generate random double [-0.5, 0.5]
				double random_num_x = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
				double random_num_y = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
				double random_num_z = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );

				// jitter sample with random numbers
				double gridPoxX = min_bounds[VX] + x * dx + random_num_x * dx;
				double gridPosY = min_bounds[VY] + y * dx + random_num_y * dx;
				double gridPosZ = min_bounds[VZ] + z * dx + random_num_z * dx;

				// add position to source_pos_list
				source_pos_list.push_back( vec3( gridPoxX,
												 gridPosY,
												 gridPosZ ) );
			}
		}





		// TODO: melting logic
		// TODO: should probably move this melting logic into its own method

		////////////////////////////////////////////////////
		// melt mesh
		////////////////////////////////////////////////////

		//if ( melting_rate != 0 ) {

		//	// TODO: delete the memory allocated for verticesValueArray after we're done with it

		//	double* verticesValueArray = new double[phi_grid.a.size()];
		//	//After the Marching Cube Method, the points of new mesh will be stored in marchingCubePointList
		//	vector<vec3> marchingCubePointList;
		//	//After the Marching Cube Method, the face list of new mesh will be stored in triangleList
		//	vector<vector<int>> triangleList; 
		//	vector<vector<int>> temp;
		//	//The origin position is the left bottom corner of the box
		//	vec3 originPos(min_box[0], min_box[1], min_box[2]);

		//	//Compute the melting thickness, the higher value of meltingThickness means the higher speed of melting
		//	double meltingThickness = maxSideLength * meltingSpeed / 100.0f;

		//	//Convert the sign of each grid's distance function
		//	for(int i = 0 ;i < phi_grid.a.size(); i++)
		//	{
		//		verticesValueArray[i] = -phi_grid.a[i];
		//	}
		//	/////////////////////////////////////////////////////////////////////////////////
		//	//From returning distance function with marching cube method to create new mesh//
		//	/////////////////////////////////////////////////////////////////////////////////
		//	for(int z = 0 ; z < sizes[2]-1; z++)
		//	for(int y = 0 ; y < sizes[1]-1; y++)
		//	for(int x = 0 ; x < sizes[0]-1; x++)
		//	{
		//		marchingCube(meltingThickness, x, y, z,//meltingThickness
		//						sizes[0]-1, sizes[1]-1, sizes[2]-1, dx,
		//						verticesValueArray, originPos, 
		//						&marchingCubePointList,
		//						&triangleList,
		//						&temp);
		//	}
		//		
		//	////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//	//Use marchingCubePointList & triangleList to create the OBJ file and import this OBJ file as a new mesh
		//	//After importing the OBJ file, it should delete the old mesh, but it will crash MAYA with some weird reason
		//	//TODO fix the problem
		//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//	stringstream ss;
		//	string tmp;
		//	string newMeshNameStr = "meltingmesh";

		//	ss<<objectCount;
		//	ss>>tmp;
		//	newMeshNameStr += tmp;
		//	ss.clear();

		//	char* newMeshNameChar=(char *)newMeshNameStr.c_str();
		//	MString newMeshNameMStr = newMeshNameChar;
		//	string test2 = newMeshNameMStr.asChar();
		//	//MString fileName = mllPath + "/meltingmesh" +(int)objectCount + ".obj";
		//	MString fileName = mllPath + "/" + newMeshNameMStr + ".obj";
		//	string fileNameStr = fileName.asChar();
		//	createObjFile(fileNameStr, marchingCubePointList, triangleList);
		//	MString importObj = "file -import -type OBJ -ra true -namespace \"" + newMeshNameMStr + "\" -options \"mo=1\"  -pr -loadReferenceDepth \"all\" \"" + fileName + "\"";
		//	string test = importObj.asChar();
		//	MGlobal::executeCommand(importObj);
		//}





	}

	// copy source_pos_list into the m_source_pos_list member variable
	m_source_pos_list.clear();
	m_source_pos_list = source_pos_list;
}