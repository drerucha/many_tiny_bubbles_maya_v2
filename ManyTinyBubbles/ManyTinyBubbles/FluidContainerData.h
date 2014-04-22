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

	// fraction field methods
	void resetFractionField( const float& default_val );
	unsigned int convert3dIndexToLinearIndex( const unsigned int& x, const unsigned int& y, const unsigned int& z ) const;
	void setFractionFieldAtXYZ( const float& val, const unsigned int& x, const unsigned int& y, const unsigned int& z );
	float getFractionFieldAtXYZ( const unsigned int& x, const unsigned int& y, const unsigned int& z ) const;

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

	float *m_fraction_field_list;

	int m_res_x;
	int m_res_y;
	int m_res_z;

	unsigned int num_voxels;

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