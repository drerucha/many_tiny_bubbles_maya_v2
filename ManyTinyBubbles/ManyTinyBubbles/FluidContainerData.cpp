#include "FluidContainerData.h"

#include <maya/MGlobal.h>
#include <maya/MItSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>

#include "Convenience.h"

#define _USE_MATH_DEFINES
#include <math.h>


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

FluidContainerData::FluidContainerData()
{
	m_initial_voxel_densities_have_been_stored = false;
}

FluidContainerData::~FluidContainerData()
{
	delete[] m_fraction_field_list;
	delete[] m_density_list;
	delete[] m_distance_function;
}


////////////////////////////////////////////////////
// initialize member variables
////////////////////////////////////////////////////
void FluidContainerData::init( MString fluid_container_name )
{
	m_name = fluid_container_name;

	MIntArray fluid_container_res_array = Convenience::getAttributeIntArray( fluid_container_name, MString( "resolution" ) );
	MDoubleArray fluid_container_dim_array = Convenience::getAttributeDoubleArray( fluid_container_name, MString( "dimensions" ) );

	// get fluid shape parent to retrieve translation attributes of fluid container
	MString fluid = Convenience::getParent( fluid_container_name );
	MDoubleArray fluid_container_translation_array = Convenience::getAttributeDoubleArray( fluid, MString( "translate" ) );

	m_res_x = fluid_container_res_array[VX];
	m_res_y = fluid_container_res_array[VY];
	m_res_z = fluid_container_res_array[VZ];

	num_voxels = m_res_x * m_res_y * m_res_z;

	m_dim_x = fluid_container_dim_array[VX];
	m_dim_y = fluid_container_dim_array[VY];
	m_dim_z = fluid_container_dim_array[VZ];

	m_trans_x = fluid_container_translation_array[VX];
	m_trans_y = fluid_container_translation_array[VY];
	m_trans_z = fluid_container_translation_array[VZ];

	m_cell_size_x = m_dim_x / m_res_x;
	m_cell_size_y = m_dim_y / m_res_y;
	m_cell_size_z = m_dim_z / m_res_z;

	// initialize fraction field data structure
	m_fraction_field_list = new double[ num_voxels ];

	// initialize constant used in computing bubble volumes
	m_sphere_volume_constant = 4.0 / 3.0 * M_PI;

	m_density_list = new double[ num_voxels ];
	m_distance_function = new double[ num_voxels ];

	storeCurrentVoxelDensities();
}


////////////////////////////////////////////////////
// get velocity field from Maya fluid container with name m_name
////////////////////////////////////////////////////
void FluidContainerData::updateVelocityField()
{
	MDoubleArray attr;
	MString cmd = "getFluidAttr -attribute \"velocity\" " + m_name;
	MGlobal::executeCommand( cmd, attr );

	m_velocity_field = attr;
}


////////////////////////////////////////////////////
// pos / array index conversion methods
////////////////////////////////////////////////////

unsigned int FluidContainerData::convertPosToLinearIndex( const vec3& pos ) const
{
	unsigned int voxel_index_x, voxel_index_y, voxel_index_z;

	convertWorldPosToGridIndices( pos,
								  voxel_index_x,
								  voxel_index_y,
								  voxel_index_z );

	return convert3dIndexToLinearIndex( voxel_index_x,
										voxel_index_y,
										voxel_index_z );
}

void FluidContainerData::convertWorldPosToGridIndices( const vec3&		pos,
													   unsigned int&	index_x,
													   unsigned int&	index_y,
													   unsigned int&	index_z ) const
{
	index_x = ( int )( ( pos[VX] + ( m_dim_x / 2.0f ) - m_trans_x ) / m_cell_size_x );
	index_y = ( int )( ( pos[VY] + ( m_dim_y / 2.0f ) - m_trans_y ) / m_cell_size_y );
	index_z = ( int )( ( pos[VZ] + ( m_dim_z / 2.0f ) - m_trans_z ) / m_cell_size_z );
}

unsigned int FluidContainerData::convert3dIndexToLinearIndex( const unsigned int& x,
															  const unsigned int& y,
															  const unsigned int& z ) const
{
	// col + row + stack
	return ( x ) + ( y * m_res_x ) + ( z * m_res_x * m_res_y );
}


////////////////////////////////////////////////////
// get fluid velocity of cell at passed-in index
////////////////////////////////////////////////////
vec3 FluidContainerData::getVelocityOfVoxelAtPos( const vec3& pos ) const
{
	// TODO: return velocity at pos; currently, this method returns the fluid velocity of the cell pos lies inside
	// TODO: bounds checking, return 0 or something if pos is outside fluid container
	// TODO: make sure this logic is correct for indexing grid faces b/c velocity is stored at cell faces, not at cell centers

	unsigned int vec3_index = convertPosToLinearIndex( pos );

	// get velocity of cell that matches the computed indices
	if ( vec3_index > 0 && vec3_index < num_voxels ) {
		double vel_x, vel_y, vel_z;
		vel_x = m_velocity_field[3 * vec3_index + 0]; // add 0 for x component
		vel_y = m_velocity_field[3 * vec3_index + 1]; // add 1 for y component
		vel_z = m_velocity_field[3 * vec3_index + 2]; // add 2 for z component

		return vec3( vel_x, vel_y, vel_z );
	}
	else {
		Convenience::printInScriptEditor( MString( "ERROR: velocity field index out of bounds in FluidContainerData::getVelocityOfVoxelAtPos" ) );
		return vec3( 0.0, 0.0, 0.0 );
	}
}


////////////////////////////////////////////////////
// fraction field methods
////////////////////////////////////////////////////

void FluidContainerData::resetFractionField()
{
	// set every element in m_fraction_field_list to 1.0f
	for ( unsigned int i = 0; i < num_voxels; ++i ) {
		m_fraction_field_list[i] = 1.0f;
	}
}

void FluidContainerData::reduceFractionFieldOfVoxelAtPos( const vec3&	bubble_pos,
														  const double&	bubble_radius )
{
	// TODO: adjust so non-square voxels work

	unsigned int linear_index_of_voxel = convertPosToLinearIndex( bubble_pos );

	// reduce fraction field of the voxel the bubble with radius bubble_radius at bubble_pos is inside
	if ( linear_index_of_voxel > 0 && linear_index_of_voxel < num_voxels ) {
		m_fraction_field_list[linear_index_of_voxel] -= m_sphere_volume_constant * pow( bubble_radius / m_cell_size_x, 3 );

		if ( m_fraction_field_list[linear_index_of_voxel] < 0.0 ) {
			m_fraction_field_list[linear_index_of_voxel] = 0.0;
		}
	}
	else {
		Convenience::printInScriptEditor( MString( "ERROR: fraction field index out of bounds in FluidContainerData::reduceFractionFieldOfVoxelAtPos" ) );
	}
}

double FluidContainerData::getFractionFieldOfVoxelAtPos( const vec3& pos ) const
{
	unsigned int index = convertPosToLinearIndex( pos );

	if ( index > 0 && index < num_voxels ) {
		return m_fraction_field_list[index];
	}
	else {
		Convenience::printInScriptEditor( MString( "ERROR: fraction field index out of bounds in FluidContainerData::getFractionFieldAtXYZ" ) );
		return 0.0f;
	}
}

//double FluidContainerData::getFractionFieldAtXYZ( const unsigned int& x, const unsigned int& y, const unsigned int& z ) const
//{
//	unsigned int index = convert3dIndexToLinearIndex( x, y, z );
//
//	if ( index > 0 && index < num_voxels ) {
//		return m_fraction_field_list[index];
//	}
//	else {
//		Convenience::printInScriptEditor( MString( "ERROR: fraction field index out of bounds in FluidContainerData::getFractionFieldAtXYZ" ) );
//		return 0.0f;
//	}
//}

//void FluidContainerData::setFractionFieldAtXYZ( const float& val, const unsigned int& x, const unsigned int& y, const unsigned int& z )
//{
//	unsigned int index = convert3dIndexToLinearIndex( x, y, z );
//
//	if ( index > 0 && index < num_voxels ) {
//		m_fraction_field_list[index] = val;
//	}
//	else {
//		Convenience::printInScriptEditor( MString( "ERROR: fraction field index out of bounds in FluidContainerData::setFractionFieldAtXYZ" ) );
//	}
//}


////////////////////////////////////////////////////
// check if pos is outside bounds of fluid container
////////////////////////////////////////////////////
bool FluidContainerData::posIsOutsideFluidContainer( const vec3& pos ) const
{
	if ( pos[VX] - m_trans_x < -m_dim_x / 2.0f ||
		 pos[VX] - m_trans_x >  m_dim_x / 2.0f ||
		 pos[VY] - m_trans_y < -m_dim_y / 2.0f ||
		 pos[VY] - m_trans_y >  m_dim_y / 2.0f ||
		 pos[VZ] - m_trans_z < -m_dim_z / 2.0f ||
		 pos[VZ] - m_trans_z >  m_dim_z / 2.0f )
	{
		return true;
	}
	else {
		return false;
	}
}


////////////////////////////////////////////////////
// check if pos is outside bounds of fluid by density value
////////////////////////////////////////////////////
bool FluidContainerData::isPosOutsideFluidViaDensity( const vec3& pos ) const
{

	if ( pos[VX] - m_trans_x < -m_dim_x / 2.0 ||
		 pos[VX] - m_trans_x >  m_dim_x / 2.0 ||
		 pos[VZ] - m_trans_z < -m_dim_z / 2.0 ||
		 pos[VZ] - m_trans_z >  m_dim_z / 2.0 )
	{
		return true;
	}
	else if ( getDensityOfVoxelAtPos( pos ) < 0.9 ) {
		unsigned int index_y = ( int )( ( pos[VY] + ( m_dim_y / 2.0f ) - m_trans_y ) / m_cell_size_y );
		double distance_grid_bottom = pos[VY] + ( m_dim_y / 2.0f ) - m_trans_y - index_y * m_cell_size_y;
		double distance_grid_bottom_percentage = distance_grid_bottom / m_cell_size_y;

		if ( distance_grid_bottom_percentage > getDensityOfVoxelAtPos( pos ) ) {
			return true;
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}
}


////////////////////////////////////////////////////
// get fluid density of cell at passed-in index
////////////////////////////////////////////////////
double FluidContainerData::getDensityOfVoxelAtPos( const vec3& pos ) const
{
	unsigned int vec3_index = convertPosToLinearIndex( pos );

	// get density of cell that matches the computed indices
	if ( vec3_index > 0 && vec3_index < num_voxels ) {
		return m_density_list[vec3_index];
	}
	else {
		Convenience::printInScriptEditor( MString( "ERROR: density index out of bounds in FluidContainerData::getDensityOfVoxelAtPos" ) );
		return 0.0;
	}
}


////////////////////////////////////////////////////
// check if pos is outside bounds of fluid by level set signed distance function
////////////////////////////////////////////////////
bool FluidContainerData::isPosOutsideFluidViaLevelSet( const vec3& pos ) const
{
	Vec3f particle_position( ( float )pos[0], ( float )pos[1], ( float )pos[2] );

	Array3f signed_distance_values;
	make_level_set3( m_fluid_mesh_face_list,
					 m_fluid_mesh_vert_list,
					 particle_position,
					 1,
					 1, 1, 1,
					 signed_distance_values );

	if ( signed_distance_values.a[0] > 0.0 ) {
		return true;
	}
	else {
		return false;
	}
	return true;
}


////////////////////////////////////////////////////
// posIsUnderBoundary()
////////////////////////////////////////////////////
int FluidContainerData::posIsUnderBoundary( const vec3& pos ) const
{
	if ( getDistanceFunctionOfVoxelAtPos( pos ) > m_cell_size_x ) {
		return 1;
	}
	else if ( getDistanceFunctionOfVoxelAtPos( pos ) <= m_cell_size_x &&
			  getDistanceFunctionOfVoxelAtPos( pos ) >= -m_cell_size_x )
	{
			return 0;
	}
	else {
		return -1;
	}
}


////////////////////////////////////////////////////
// get distance function of cell
////////////////////////////////////////////////////
double FluidContainerData::getDistanceFunctionOfVoxelAtPos( const vec3& pos ) const
{
	unsigned int vec3_index = convertPosToLinearIndex( pos );

	// get distance function of cell that matches the computed indices
	if ( vec3_index > 0 && vec3_index < num_voxels ) {
		return m_distance_function[vec3_index];
	}
	else {
		Convenience::printInScriptEditor( MString( "ERROR: distance function index out of bounds in FluidContainerData::getDistanceFunctionOfVoxelAtPos" ) );
		return 0.0;
	}
}


////////////////////////////////////////////////////
// updateDensityField()
////////////////////////////////////////////////////
void FluidContainerData::updateDensityField()
{
	MDoubleArray attr;
	MString cmd = "getFluidAttr -attribute \"density\" " + m_name;
	MGlobal::executeCommand( cmd, attr );

	for ( unsigned int i = 0; i < num_voxels; ++i ) {
		m_density_list[i] = attr[i];
	}
}


////////////////////////////////////////////////////
// prepareFluidMeshForLevelSetMethod()
////////////////////////////////////////////////////
void FluidContainerData::prepareFluidMeshForLevelSetMethod( MString fluid_polygon_name )
{
	////////////////////////////////////////////////////
	// first, get Maya mesh object from mesh name
	////////////////////////////////////////////////////

	// select Maya object by name
	MSelectionList maya_sel_list;
	MGlobal::getSelectionListByName( fluid_polygon_name, maya_sel_list );

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

	
		////////////////////////////////////////////////////
		// create face list
		// use Vec3ui for compatibility with SDFGen level set library
		////////////////////////////////////////////////////
		m_fluid_mesh_face_list.clear();
		for ( unsigned int i = 0; i < triangle_vertex_num; i += 3 ) {
			m_fluid_mesh_face_list.push_back( Vec3ui( triangle_vertex_indices[i],
													  triangle_vertex_indices[i + 1],
													  triangle_vertex_indices[i + 2] ) );
		}


		////////////////////////////////////////////////////
		// create vertex list
		// use Vec3f for compatibility with SDFGen level set library
		////////////////////////////////////////////////////
		MFloatPointArray mesh_vertices;
		mesh_surface.getPoints(mesh_vertices, MSpace::kWorld);
		m_fluid_mesh_vert_list.clear();
		for ( unsigned int i = 0; i < mesh_vertices.length(); ++i ) {
			MFloatPoint vertex = mesh_vertices[i];
			Vec3f new_vert( vertex[VX],
							vertex[VY],
							vertex[VZ] );
			m_fluid_mesh_vert_list.push_back( new_vert );
		}
	}
}


////////////////////////////////////////////////////
// update the current Maya fluid level set signed distance function
////////////////////////////////////////////////////
void FluidContainerData::updateLevelSetSignedDistanceFunction( MString fluid_polygon_name ) const
{
	////////////////////////////////////////////////////
	// first, get Maya mesh object from mesh name
	////////////////////////////////////////////////////

	// select Maya object by name
	MSelectionList maya_sel_list;
	MGlobal::getSelectionListByName( fluid_polygon_name, maya_sel_list );

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

		double bottom_left_corner_x = m_trans_x - ( m_dim_x / m_res_x ) * ( m_res_x / 2.0 - 0.5 );
		double bottom_left_corner_y = m_trans_y - ( m_dim_y / m_res_y ) * ( m_res_y / 2.0 - 0.5 );
		double bottom_left_corner_z = m_trans_z - ( m_dim_z / m_res_z ) * ( m_res_z / 2.0 - 0.5 );

		Vec3f leftBottomCorner( ( float )bottom_left_corner_x,
								( float )bottom_left_corner_y,
								( float )bottom_left_corner_z );

		float dx = ( float )m_cell_size_x;
		min_bounds += Vec3f( ( float )m_cell_size_x,
							 ( float )m_cell_size_y,
							 ( float )m_cell_size_z );


		////////////////////////////////////////////////////
		// compute signed distance function
		////////////////////////////////////////////////////

		// number of sample points per x, y, z direction?
		Vec3ui sizes = Vec3ui( m_res_x, m_res_y, m_res_z );

		// Array3f is a data structure defined in array3.h of SDFGen library
		Array3f signed_distance_values;

		// make_level_set3 is a method in makelevelset3.cpp of SDFGen library
		make_level_set3( face_list,
						 vert_list,
						 leftBottomCorner,
						 dx,
						 sizes[VX], sizes[VY], sizes[VZ],
						 signed_distance_values );

		for ( unsigned int i = 0; i < num_voxels; ++i ) {
			m_distance_function[i] = signed_distance_values.a[i];
		}
	}
}


////////////////////////////////////////////////////
// populate m_initial_voxel_densities with current voxel densities of fluid
// should only be called one time when fluid is initialized
////////////////////////////////////////////////////
void FluidContainerData::storeCurrentVoxelDensities()
{
	if ( !m_initial_voxel_densities_have_been_stored ) {
		m_initial_voxel_densities_have_been_stored = true;

		MDoubleArray attr;
		MString cmd = "getFluidAttr -attribute \"density\" " + m_name;
		MGlobal::executeCommand( cmd, attr );

		m_voxel_densities = attr;
	}
}


////////////////////////////////////////////////////
// fill fluid container voxels with densities stored in m_voxel_densities
////////////////////////////////////////////////////
void FluidContainerData::resetDensity()
{
	for ( unsigned int index_y = 0; index_y < m_res_y; ++index_y ) {
		for ( unsigned int index_z = 0; index_z < m_res_z; ++index_z ) {
			for ( unsigned int index_x = 0; index_x < m_res_x; ++index_x ) {

				int index = convert3dIndexToLinearIndex( index_x, index_y, index_z );
				MString mstr_density_val = Convenience::convertDoubleToMString( m_voxel_densities[index] );

				MString mstr_index_x = Convenience::convertIntToMString( index_x );
				MString mstr_index_y = Convenience::convertIntToMString( index_y );
				MString mstr_index_z = Convenience::convertIntToMString( index_z );

				//setFluidAttr -at "density"
				//			 -fv 1.0
				//			 -xIndex $index_x
				//			 -yIndex $index_y
				//			 -zIndex $index_z
				//			 $fluid_container;

				MString cmd = "setFluidAttr -attribute \"density\" -floatValue " + mstr_density_val + " -xIndex " + mstr_index_x + " -yIndex " + mstr_index_y + " -zIndex " + mstr_index_z + " " + m_name;
				MGlobal::executeCommand( cmd );
			}
		}
	}
}