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

 
class ManyTinyBubbles : public MPxNode
{
public:
						ManyTinyBubbles();
	virtual				~ManyTinyBubbles();

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static  void*		creator();
	static  MStatus		initialize();

public:

	// There needs to be a MObject handle declared for each attribute that
	// the node will have.  These handles are needed for getting and setting
	// the values later.

	static MObject input;	// test
	static MObject output;	// test

	static MObject time;
	static MObject emitter_mesh_name;
	static MObject fluid_container_name;
	static MObject emission_rate;
	static MObject scattering_frequency;
	static MObject scattering_coefficient;
	static MObject breakup_frequency;
	static MObject bubble_size_min;
	static MObject bubble_size_max;

	// The typeid is a unique 32bit indentifier that describes this node.
	// It is used to save and retrieve nodes of this type from the binary
	// file format.  If it is not unique, it will cause file IO problems.
	static	MTypeId		id;
};

#endif
