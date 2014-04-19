#pragma once

#include <maya/MString.h>
#include <maya/MDoubleArray.h>

#include "vec.h"


class FluidContainerData
{

////////////////////////////////////////////////////
// methods
////////////////////////////////////////////////////

public:

	FluidContainerData( void );
	~FluidContainerData( void );

	void init( MString fluid_container_name );

	void updateVelocityField( void );

	vec3 getVelocityAtPos( const vec3& pos ) const;

private:

	void convertWorldPosToGridIndices( const vec3&		pos,
									   unsigned int&	index_x,
									   unsigned int&	index_y,
									   unsigned int&	index_z ) const;


////////////////////////////////////////////////////
// members
////////////////////////////////////////////////////

private:

	MString m_name;

	MDoubleArray m_velocity_field;

	int m_res_x;
	int m_res_y;
	int m_res_z;

	double m_dim_x;
	double m_dim_y;
	double m_dim_z;

	double m_trans_x;
	double m_trans_y;
	double m_trans_z;

	double m_cell_size_x;
	double m_cell_size_y;
	double m_cell_size_z;

};