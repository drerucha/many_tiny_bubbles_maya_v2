#pragma once

#include <maya/MString.h>
#include <maya/MDoubleArray.h>

#include "vec.h"
#include "SDFGen/makelevelset3.h"


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

	vec3 getVelocityOfVoxelAtPos( const vec3& pos ) const;
	//vec3 getVelocityOfVoxel( const int& x, const int& y, const int& z ) const;

	// fraction field methods
	void resetFractionField();
	void reduceFractionFieldOfVoxelAtPos( const vec3&	bubble_pos,
										  const double&	bubble_radius );
	double getFractionFieldOfVoxelAtPos( const vec3& pos ) const;
	//double getFractionFieldAtXYZ( const unsigned int& x, const unsigned int& y, const unsigned int& z ) const;
	//void setFractionFieldAtXYZ( const float& val, const unsigned int& x, const unsigned int& y, const unsigned int& z );

	bool posIsOutsideFluidContainer( const vec3& pos ) const;

	void prepareFluidMeshForLevelSetMethod( MString fluid_polygon_name );

	//bool posIsOutsideFluidLevel_Density( const vec3& pos ) const;
	//bool posIsOutsideFluidLevel_LevelSet( const vec3& pos ) const;
	//int posIsUnderBoundary( const vec3& pos ) const;
	//void updateFluidLevelSetDistanceFunction( MString fluid_polygon_name ) const;
	//void updateDensityField();
	//double getDensityOfVoxelAtPos( const vec3& pos ) const;
	//double getDistanceFunctionOfVoxelAtPos( const vec3& pos ) const;

private:

	// pos / array index conversion methods
	unsigned int convertPosToLinearIndex( const vec3& pos ) const;
	void convertWorldPosToGridIndices( const vec3&		pos,
									   unsigned int&	index_x,
									   unsigned int&	index_y,
									   unsigned int&	index_z ) const;
	unsigned int convert3dIndexToLinearIndex( const unsigned int& x,
											  const unsigned int& y,
											  const unsigned int& z ) const;


////////////////////////////////////////////////////
// members
////////////////////////////////////////////////////

private:

	MString m_name;

	MDoubleArray m_velocity_field;

	// how much water is in each voxel
	double *m_fraction_field_list;

	double *m_density_list;
	double *m_distance_function;

	std::vector<Vec3ui> m_fluid_mesh_face_list;
	std::vector<Vec3f> m_fluid_mesh_vert_list;

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

	double m_sphere_volume_constant;
};