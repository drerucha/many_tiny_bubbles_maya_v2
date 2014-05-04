#include "EmitterData.h"

#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>

#include "SDFGen/makelevelset3.h"
#include "Convenience.h"
#include "GlobalState.h"

// TODO: check if mesh is a sphere, and generate emission positions as appropriate


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
void EmitterData::init( int emit_rate )
{
	m_emission_rate = emit_rate;
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

	// get all selected mesh names
	std::vector<std::string> selectedObjectNames = GlobalState::getSelectedObject();
	for ( int objectCount = 0; objectCount< selectedObjectNames.size(); ++objectCount ) {

		// retrieve the selected mesh object
		MString strSelectObjName = ( char* )selectedObjectNames[objectCount].c_str();

		// if this object does not exist, delete the name from list
		int exist;
		MGlobal::executeCommand( "objExists " + strSelectObjName, exist );
		if ( exist == 0 ) {
			GlobalState::deleteSelectedObject( Convenience::convertMStringToStdString(strSelectObjName) );
			continue;
		}

		MStringArray emitterMeshInfoArray;

		// script to get the material name of this mesh
		MGlobal::executeCommand( "listSets -t 1 -ets -o " + strSelectObjName, emitterMeshInfoArray );
		MString emitterMeshMaterial = emitterMeshInfoArray[0];
			

		////////////////////////////////////////////////////
		// first, get Maya mesh object from mesh name
		////////////////////////////////////////////////////

		// select Maya object by name
		MSelectionList maya_sel_list;
		MGlobal::getSelectionListByName( strSelectObjName, maya_sel_list );

		// get path to a mesh DAG node
		MDagPath dag_node_path;
		maya_sel_list.getDagPath( 0, dag_node_path );

		// if dag_node_path points to a transform node instead of a shape node, try to extend the DAG path to reach a shape node
		bool node_has_shape = true;
		if ( dag_node_path.apiType() == MFn::kTransform ) {
			MStatus stat = dag_node_path.extendToShape();
			if ( stat != MStatus::kSuccess ) {
				node_has_shape = false;
				GlobalState::deleteSelectedObject( Convenience::convertMStringToStdString(strSelectObjName) );
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

			// compute reasonable dx for computing level sets
			// dx represents the level set resolution
			float bounding_box_volume = bounding_box_size_x * bounding_box_size_y * bounding_box_size_z;
			float voxel_volume = bounding_box_volume / voxel_num; // voxel_num is a method argument
			//float dx = pow( voxel_volume, ( 1.0f / 3.0f ) );
			m_level_set_dx = pow( voxel_volume, ( 1.0f / 3.0f ) );

			// add padding around bounding box
			Vec3f unit( 1.0f, 1.0f, 1.0f );
			min_bounds -= 0.001f * unit;
			max_bounds += 0.001f * unit;


			////////////////////////////////////////////////////
			// compute signed distance function
			////////////////////////////////////////////////////

			// number of sample points per x, y, z direction?
			Vec3ui sizes = Vec3ui( ( max_bounds - min_bounds ) / ( float )m_level_set_dx ) + Vec3ui( 2, 2, 2 );

			// Array3f is a data structure defined in array3.h of SDFGen library
			Array3f signed_distance_values;

			// make_level_set3 is a method in makelevelset3.cpp of SDFGen library
			make_level_set3( face_list,
							 vert_list,
							 min_bounds,
							 ( float )m_level_set_dx,
							 sizes[VX],
							 sizes[VY],
							 sizes[VZ],
							 signed_distance_values );

			// create bubble generation positions from signed distance function comuted in make_level_set3
			// signed_distance_values.a is an ArrayT
			for ( unsigned int i = 0; i < signed_distance_values.a.size(); ++i ) {

				// bubble generation position will be on outer surface of mesh or within some threshold above it
				if ( signed_distance_values.a[i] >= 0.0f && signed_distance_values.a[i] < m_level_set_dx / 2.0f ) {

					// TODO: ask what ni, nj, nk are for
					// TODO: ask about this position computation logic

					// convert linear signed distance array index into 3D indices
					int z = i / signed_distance_values.ni / signed_distance_values.nj;
					int y = ( i - z * signed_distance_values.ni * signed_distance_values.nj ) / signed_distance_values.ni;
					int x = i - z * signed_distance_values.ni * signed_distance_values.nj - y * signed_distance_values.ni;

					double gridPoxX = min_bounds[VX] + x * m_level_set_dx;
					double gridPosY = min_bounds[VY] + y * m_level_set_dx;
					double gridPosZ = min_bounds[VZ] + z * m_level_set_dx;

					// add position to source_pos_list
					source_pos_list.push_back( vec3( gridPoxX,
													 gridPosY,
													 gridPosZ ) );
				}
			}


			// TODO: melting logic
			// TODO: should probably move this melting logic into its own method
		}
	}

	// copy source_pos_list into the m_source_pos_list member variable
	m_source_pos_list.clear();
	m_source_pos_list = source_pos_list;
}


////////////////////////////////////////////////////
// generate bubbles on mesh or sphere surface
// return new bubble positions, velocities, and radius group
////////////////////////////////////////////////////
void EmitterData::generateBubbles( const std::vector<double>&	bubble_radii_list,
								   std::vector<vec3>&			new_bubble_positions,
								   std::vector<vec3>&			new_bubble_velocities,
								   std::vector<unsigned int>&	new_bubble_radius_group ) const
{
	// TODO: make bubble generate rate dependent on user-defined emission rate

	int total_bubbles_per_frame = m_emission_rate;
	const vec3 INITIAL_BUBBLE_VELOCITY( 0.0, 1.0, 0.0 );

	// if emitter is a mesh
	if ( m_source_pos_list.size() != 0 ) {

		// iterate through bubble radius groups
		for ( unsigned int k = 0; k < bubble_radii_list.size(); ++k ) {

			// divide the total bubbles per frame with frame number, and round the result
			unsigned int bubbles_in_this_radii = ( int )std::floor( ( float )total_bubbles_per_frame / ( ( float )bubble_radii_list.size() - k ) + 0.5f );
			total_bubbles_per_frame -= bubbles_in_this_radii;

			if ( bubbles_in_this_radii == 0 ) {
				continue;
			}

			// if the bubbles in this radius group are more than the source postions
			if ( bubbles_in_this_radii > m_source_pos_list.size() ) {

				// iterate through every available emission position
				for ( unsigned int i = 0; i < m_source_pos_list.size(); ++i ) {

					// compute the bubbles in each position and round the number
					unsigned int bubbles_per_position = ( int )std::floor( ( float )bubbles_in_this_radii / ( float )m_source_pos_list.size() + 0.5f ) ;
					bubbles_in_this_radii -= bubbles_per_position;

					// each position has multiple bubbles, so iterate to generate bubbles
					for ( unsigned int p = 0; p < bubbles_per_position; ++p ) {
						vec3 bubble_position = m_source_pos_list[i];

						double random_num_x = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
						double random_num_y = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
						double random_num_z = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );

						// jitter the position
						bubble_position += vec3( m_level_set_dx * random_num_x,
												 m_level_set_dx * random_num_y,
												 m_level_set_dx * random_num_z );

						new_bubble_positions.push_back( bubble_position );
						new_bubble_velocities.push_back( INITIAL_BUBBLE_VELOCITY );
						new_bubble_radius_group.push_back( k );
					}
				}

			}
			// if the bubbles in this radius group are less than the source postions
			else {

				// iterate through every available emission position
				// but because the postions are more than the bubbles, bubbles are generated on period positions
				for ( unsigned int i = 0; i < bubbles_in_this_radii; ++i ) {

					//get random index from the total position list
					int position_index = Convenience::generateRandomIntInclusive( 0, ( int )m_source_pos_list.size() - 1 );
					vec3 bubble_position = m_source_pos_list[position_index];

					double random_num_x = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
					double random_num_y = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
					double random_num_z = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );

					// jitter the position
					bubble_position += vec3( m_level_set_dx * random_num_x,
											 m_level_set_dx * random_num_y,
											 m_level_set_dx * random_num_z );

					new_bubble_positions.push_back( bubble_position );
					new_bubble_velocities.push_back( INITIAL_BUBBLE_VELOCITY );
					new_bubble_radius_group.push_back( k );
				}
			}
		}
	}
}