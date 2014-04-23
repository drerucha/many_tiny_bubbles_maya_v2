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

	void setSourcePosList( const unsigned int& cubicNum,
						   const unsigned int& meltingSpeed );

private:

	MString m_name;

	// list of positions on mesh where bubbles can emit from
	std::vector<vec3> m_source_pos_list;

};

// TODO: move emitter attributes from BubbleData to here
// TODO: make mesh emitter name a member variable of this class
// TODO: create method to generate generation positions from mesh