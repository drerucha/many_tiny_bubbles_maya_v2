#pragma once

#include <maya/MString.h>

#include "vec.h"
#include <vector>


class EmitterData
{
public:
	EmitterData( void );
	~EmitterData( void );

	void init( MString name );

	void createEmissionPositionsOnMesh( const unsigned int& voxel_num,
										const unsigned int& melting_rate );

	void generateBubbles( const std::vector<double>&	bubble_radii_list,
						  std::vector<vec3>&			new_bubble_positions,
						  std::vector<vec3>&			new_bubble_velocities,
						  std::vector<unsigned int>&	new_bubble_radius_group ) const;

private:

	MString m_name;

	// list of positions on mesh where bubbles can emit from
	std::vector<vec3> m_source_pos_list;

	vec3 m_sphere_center;
	double m_sphere_radius;

};

// TODO: move emitter attributes from BubbleData to here
// TODO: create method to generate generation positions from mesh
// TODO: add the bounding box as a member variable?