#ifndef _ManyTinyBubblesNode
#define _ManyTinyBubblesNode
//
// Copyright (C) CIS660
// 
// File: ManyTinyBubblesNode.h
//
// Dependency Graph Node: ManyTinyBubbles
//
// Author: Maya Plug-in Wizard 2.0
//

#include <maya/MPxNode.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MTypeId.h>

#include "BubbleData.h"
#include "FluidContainerData.h"

 
class ManyTinyBubbles : public MPxNode
{

////////////////////////////////////////////////////
// methods
////////////////////////////////////////////////////

public:

	ManyTinyBubbles();
	virtual ~ManyTinyBubbles();

	virtual MStatus	compute( const MPlug& plug, MDataBlock& data );

	static void* creator();
	static MStatus initialize();

private:

// for debugging
void testCode( const MString& str ) const;

MStatus	createBubbles( const MTime& time,
					   const float& step_size );

void computeFractionField( void );

void advectParticles( const float& dt );

void reset();


////////////////////////////////////////////////////
// members
////////////////////////////////////////////////////

public:

	// there must be an MObject handle declared for each attribute that the node will have

	//static MObject input;
	static MObject m_output;
	static MObject m_time;
	//static MObject m_emitter_mesh;
	static MObject m_emitter_mesh_name;
	static MObject m_fluid_container_name;
	static MObject m_emission_rate;
	static MObject m_scattering_frequency;
	static MObject m_scattering_coefficient;
	static MObject m_breakup_frequency;
	static MObject m_bubble_size_min;
	static MObject m_bubble_size_max;
	static MObject m_step_size;

	// typeid is a unique 32bit indentifier that describes this node
	static	MTypeId m_id;

private:

	BubbleData m_bubbles;
	FluidContainerData m_fluid_container;

	unsigned int m_current_frame;
};

#endif
