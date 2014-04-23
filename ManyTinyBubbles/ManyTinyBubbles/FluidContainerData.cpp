#include "FluidContainerData.h"

#include <maya/MGlobal.h>

#include "Convenience.h"

#define _USE_MATH_DEFINES
#include <math.h>


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

FluidContainerData::FluidContainerData()
{
}

FluidContainerData::~FluidContainerData()
{
	delete[] m_fraction_field_list;
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

	int vec3_index = convertPosToLinearIndex( pos );

	// get velocity of cell that matches the computed indices
	double vel_x, vel_y, vel_z;
	vel_x = m_velocity_field[vec3_index + 0]; // add 0 for x component
	vel_y = m_velocity_field[vec3_index + 1]; // add 1 for y component
	vel_z = m_velocity_field[vec3_index + 2]; // add 2 for z component

	return vec3( vel_x, vel_y, vel_z );
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

	int linear_index_of_voxel = convertPosToLinearIndex( bubble_pos );

	// reduce fraction field of the voxel the bubble with radius bubble_radius at bubble_pos is inside
	m_fraction_field_list[linear_index_of_voxel] -= m_sphere_volume_constant * pow( bubble_radius / m_cell_size_x, 3 );
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