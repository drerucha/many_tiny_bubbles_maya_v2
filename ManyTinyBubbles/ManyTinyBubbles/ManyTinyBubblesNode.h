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

	// must be an MObject handle declared for each attribute that the node will have

	static MObject input;
	static MObject output;
	static MObject time;
	static MObject emitter_mesh;
	static MObject emitter_mesh_name;
	static MObject fluid_container_name;
	static MObject emission_rate;
	static MObject scattering_frequency;
	static MObject scattering_coefficient;
	static MObject breakup_frequency;
	static MObject bubble_size_min;
	static MObject bubble_size_max;

	// typeid is a unique 32bit indentifier that describes this node
	static	MTypeId id;


private:

	int fluid_container_res_x;
	int fluid_container_res_y;
	int fluid_container_res_z;

	double fluid_container_dim_x;
	double fluid_container_dim_y;
	double fluid_container_dim_z;

	double fluid_container_trans_x;
	double fluid_container_trans_y;
	double fluid_container_trans_z;

	double fluid_container_cell_size_x;
	double fluid_container_cell_size_y;
	double fluid_container_cell_size_z;

};

#endif
