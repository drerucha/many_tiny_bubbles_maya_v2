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
#include "EmitterData.h"

 
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

MStatus	createBubbles( const MTime& time, const float& step_size );

void advectParticles( const float& dt );

double computeScatteringProbabilityOfBubble( const vec3& bubble_vel, const vec3& bubble_pos ) const;

double computeAlteredAngle( void ) const;

vec3 updateBubbleVelocity( const vec3& old_vel, const double& altered_dir ) const;

void splitBubble( const vec3&	current_pos,
				  const double&	current_radius,
				  vec3&			new_pos_1,
				  vec3&			new_pos_2 ) const;


struct BubbleLocator {
	unsigned int radius_group;
	unsigned int list_index;

	BubbleLocator( unsigned int i, unsigned int j )
	{
		radius_group = i;
		list_index = j;
	}
};

struct BubbleStruct {
	vec3 pos;
	vec3 vel;
	unsigned int radius_group;

	BubbleStruct( vec3 pos,
				  vec3 vel,
				  unsigned int radius_group )
	{
		this->pos = pos;
		this->vel = vel;
		this->radius_group = radius_group;
	}
};

MStatus	simulationLoop( const MTime& time, const float& step_size ); // Danny was here
unsigned int simulationSetup( const MTime& time ); // Danny was here
void deleteEscapedBubbles( void ); // Danny was here
void deleteBubblesInList( std::vector<BubbleLocator> bubbles_to_remove ); // Danny was here
void generateMoreBubblesFromEmitter( void ); // Danny was here
void computeFractionField( void ); // Danny was here
void updateBubbleVelocities( void ); // Danny was here
void breakupBubbles( void ); // Danny was here
void reset(); // Danny was here

// TODO: add method to delete bubbles given a vector of BubbleLocator structs
// TODO: add method to add bubbles given a vector of BubbleStruct structs


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
	EmitterData m_emitter;

	unsigned int m_current_frame;
};

#endif