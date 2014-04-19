#include "FluidContainerData.h"

#include <maya/MGlobal.h>

#include "Convenience.h"


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

FluidContainerData::FluidContainerData(void)
{
}

FluidContainerData::~FluidContainerData(void)
{
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

	m_dim_x = fluid_container_dim_array[VX];
	m_dim_y = fluid_container_dim_array[VY];
	m_dim_z = fluid_container_dim_array[VZ];

	m_trans_x = fluid_container_translation_array[VX];
	m_trans_y = fluid_container_translation_array[VY];
	m_trans_z = fluid_container_translation_array[VZ];

	m_cell_size_x = m_dim_x / m_res_x;
	m_cell_size_y = m_dim_y / m_res_y;
	m_cell_size_z = m_dim_z / m_res_z;
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
// get fluid velocity of cell at passed-in index
////////////////////////////////////////////////////
vec3 FluidContainerData::getVelocityAtPos( const vec3& pos ) const
{
	// TODO: return velocity at pos; currently, this method returns the fluid velocity of the cell pos lies inside
	// TODO: bounds checking, return 0 or something if pos is outside fluid container
	// TODO: make sure this logic is correct for indexing grid faces b/c velocity is stored at cell faces, not at cell centers

	// get cell indices of fluid container for passed-in pos
	unsigned int fluid_cell_index_x, fluid_cell_index_y, fluid_cell_index_z;
	convertWorldPosToGridIndices( pos,
								  fluid_cell_index_x,
								  fluid_cell_index_y,
								  fluid_cell_index_z );

	// col + row + stack
	int vec3_index = ( fluid_cell_index_x ) + ( fluid_cell_index_y * m_res_x ) + ( fluid_cell_index_z * m_res_x * m_res_y );

	// get velocity of cell that matches the computed indices
	double vel_x, vel_y, vel_z;
	vel_x = m_velocity_field[vec3_index + 0]; // add 0 for x component
	vel_y = m_velocity_field[vec3_index + 1]; // add 1 for y component
	vel_z = m_velocity_field[vec3_index + 2]; // add 2 for z component

	return vec3( vel_x, vel_y, vel_z );
}


////////////////////////////////////////////////////
// convert world space position to indices of 
////////////////////////////////////////////////////
void FluidContainerData::convertWorldPosToGridIndices( const vec3&		pos,
													   unsigned int&	index_x,
													   unsigned int&	index_y,
													   unsigned int&	index_z ) const
{
	index_x = ( int )( ( pos[VX] + ( m_dim_x / 2.0f ) - m_trans_x ) / m_cell_size_x );
	index_y = ( int )( ( pos[VY] + ( m_dim_y / 2.0f ) - m_trans_y ) / m_cell_size_y );
	index_z = ( int )( ( pos[VZ] + ( m_dim_z / 2.0f ) - m_trans_z ) / m_cell_size_z );
}