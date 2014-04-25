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
// structs
////////////////////////////////////////////////////

private:

struct BubbleLocator {
	unsigned int	radius_group;
	unsigned int	list_index;

	BubbleLocator( unsigned int i, unsigned int j )
	{
		radius_group = i;
		list_index = j;
	}
};

struct BubbleStruct {
	vec3			pos;
	vec3			vel;
	unsigned int	radius_group;

	BubbleStruct( vec3 pos,
				  vec3 vel,
				  unsigned int radius_group )
	{
		this->pos = pos;
		this->vel = vel;
		this->radius_group = radius_group;
	}
};


////////////////////////////////////////////////////
// methods
////////////////////////////////////////////////////

public:

	// necessary for Maya node
						ManyTinyBubbles();
	virtual				~ManyTinyBubbles();
	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );
	static void*		creator();
	static MStatus		initialize();

private:

	MStatus			simulationLoop( const MTime& time, const float& step_size );
	unsigned int	simulationSetup( const MTime& time );
	void			reset();

	void			deleteEscapedBubbles( void );
	void			deleteBubblesInList( std::vector<BubbleLocator> bubbles_to_remove );

	void			generateMoreBubblesFromEmitter( void );

	void			updateBubbleVelocities( void );
	void			computeFractionField( void );
	double			computeScatteringProbabilityOfBubble( const vec3& bubble_vel, const vec3& bubble_pos ) const;
	double			computeAlteredAngle( void ) const;
	vec3			computeBubbleVelAfterScattering( const vec3& old_vel, const double& altered_dir ) const;

	void			breakupBubbles( void );
	void			splitBubble( const vec3&	current_pos,
								 const double&	current_radius,
								 vec3&			new_pos_1,
								 vec3&			new_pos_2 ) const;


////////////////////////////////////////////////////
// members
////////////////////////////////////////////////////

public:

	// there must be an MObject handle declared for each attribute that the node will have
	static MObject	m_output;
	static MObject	m_time;
	//static MObject m_emitter_mesh;
	static MObject	m_emitter_mesh_name;
	static MObject	m_fluid_container_name;
	static MObject	m_emission_rate;
	static MObject	m_scattering_frequency;
	static MObject	m_scattering_coefficient;
	static MObject	m_breakup_frequency;
	static MObject	m_bubble_size_min;
	static MObject	m_bubble_size_max;
	static MObject	m_step_size;

	// typeid is a unique 32bit indentifier that describes this node
	static MTypeId	m_id;

private:

	BubbleData			m_bubbles;
	FluidContainerData	m_fluid_container;
	EmitterData			m_emitter;

	unsigned int		m_current_frame;
};

#endif