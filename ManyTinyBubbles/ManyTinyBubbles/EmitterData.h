#pragma once

#include <maya/MString.h>
//#include <maya/MFnMesh.h>

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

private:

	MString m_name;

	//MFnMesh m_mesh;

	// list of positions on mesh where bubbles can emit from
	std::vector<vec3> m_source_pos_list;

};

// TODO: move emitter attributes from BubbleData to here
// TODO: create method to generate generation positions from mesh
// TODO: add the bounding box as a member variable?