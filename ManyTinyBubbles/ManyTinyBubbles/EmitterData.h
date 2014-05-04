//#pragma once

#ifndef _EmitterData
#define _EmitterData

#include <maya/MString.h>

#include "vec.h"
#include <vector>


class EmitterData
{
public:
	EmitterData( void );
	~EmitterData( void );

	void init( int emit_rate, double min_bubble_radius,
			   unsigned int level_set_res, unsigned int melting_rate );

	void createEmissionPositionsOnMesh( const unsigned int &voxel_num, const unsigned int &melting_rate );

	void generateBubbles( const std::vector<double>		&bubble_radii_list,
						  std::vector<vec3>				&new_bubble_positions,
						  std::vector<vec3>				&new_bubble_velocities,
						  std::vector<unsigned int>		&new_bubble_radius_group ) const;

	int marchingCube( double						&isolevel,
					  const unsigned int			&indexX,
					  const unsigned int			&indexY,
					  const unsigned int			&indexZ,
					  const unsigned int			&meshResX,
					  const unsigned int			&meshResY,
					  const unsigned int			&meshResZ,
					  double						&meshSize,
					  double						*verticesValueArray,
					  vec3							originPos,
					  std::vector<vec3>				*marchingCubePointList,
					  std::vector<std::vector<int>>	*triangleList,
					  std::vector<std::vector<int>>	*temp );

	vec3 vertexInterp( double	&isolevel,
					   vec3		p1,
					   vec3		p2,
					   double	&valueP1,
					   double	&valueP2 );

	void createObjFile( std::string fileName,
						std::vector<vec3> marchingCubePointList,
						std::vector<std::vector<int>> triangleList );

	int getSourcePosListSize( void );

	void deleteEmitterMeshesFromScene( void );
	void createObjFileFromStoredMeshData( void );

	// getters
	unsigned int getLevelSetRes( void );
	unsigned int getMeltingRate( void );

private:
	//MString m_name;

	// list of positions on mesh where bubbles can emit from
	std::vector<vec3> m_source_pos_list;

	vec3 m_sphere_center;
	double m_sphere_radius;

	double m_level_set_dx;
	int m_emission_rate;
	double m_min_bubble_radius;

	unsigned int m_level_set_res;
	unsigned int m_melting_rate;

};

#endif