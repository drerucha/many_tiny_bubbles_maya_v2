//
// Copyright 2014 (C) Danny Rerucha and Wei-Chien Tu
// 
// File: ManyTinyBubblesNode.cpp
//
// Dependency Graph Node: ManyTinyBubbles
//
// Authors: Danny Rerucha and Wei-Chien Tu
//


#include "ManyTinyBubblesNode.h"

#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

#include <maya/MGlobal.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MTime.h>


MTypeId ManyTinyBubbles::id( 0x70256 );


////////////////////////////////////////////////////
// node attributes
////////////////////////////////////////////////////

MObject ManyTinyBubbles::input;
MObject ManyTinyBubbles::output;
MObject ManyTinyBubbles::time;
MObject ManyTinyBubbles::emitter_mesh;
MObject ManyTinyBubbles::emitter_mesh_name;
MObject ManyTinyBubbles::fluid_container_name;
MObject ManyTinyBubbles::emission_rate;
MObject ManyTinyBubbles::scattering_frequency;
MObject ManyTinyBubbles::scattering_coefficient;
MObject ManyTinyBubbles::breakup_frequency;
MObject ManyTinyBubbles::bubble_size_min;
MObject ManyTinyBubbles::bubble_size_max;


////////////////////////////////////////////////////
// constructor/destructor
////////////////////////////////////////////////////

ManyTinyBubbles::ManyTinyBubbles() {}
ManyTinyBubbles::~ManyTinyBubbles() {}


////////////////////////////////////////////////////
// compute(): computes value of given output plug based on values of input attributes
//		plug - plug to compute
//		data - object that provides access to attributes for this node
////////////////////////////////////////////////////
MStatus ManyTinyBubbles::compute( const MPlug& plug, MDataBlock& data )
{
	MStatus returnStatus;


	////////////////////////////////////////////////////
	// check which output attribute we have been asked to compute
	////////////////////////////////////////////////////

	if ( plug == ManyTinyBubbles::output ) {

		////////////////////////////////////////////////////
		// get handles to input attribute we will need for computation
		////////////////////////////////////////////////////

		MDataHandle input_data						= data.inputValue( ManyTinyBubbles::input, &returnStatus );
		MDataHandle time_data						= data.inputValue( ManyTinyBubbles::time, &returnStatus );
		MDataHandle emitter_mesh_data				= data.inputValue( ManyTinyBubbles::emitter_mesh, &returnStatus );
		MDataHandle emitter_mesh_name_data			= data.inputValue( ManyTinyBubbles::emitter_mesh_name, &returnStatus );
		MDataHandle fluid_container_name_data		= data.inputValue( ManyTinyBubbles::fluid_container_name, &returnStatus );
		MDataHandle emission_rate_data				= data.inputValue( ManyTinyBubbles::emission_rate, &returnStatus );
		MDataHandle scattering_frequency_data		= data.inputValue( ManyTinyBubbles::scattering_frequency, &returnStatus );
		MDataHandle scattering_coefficient_data		= data.inputValue( ManyTinyBubbles::scattering_coefficient, &returnStatus );
		MDataHandle breakup_frequency_data			= data.inputValue( ManyTinyBubbles::breakup_frequency, &returnStatus );
		MDataHandle bubble_size_min_data			= data.inputValue( ManyTinyBubbles::bubble_size_min, &returnStatus );
		MDataHandle bubble_size_max_data			= data.inputValue( ManyTinyBubbles::bubble_size_max, &returnStatus );

		if ( returnStatus != MS::kSuccess ) {
			MGlobal::displayError( "Node ManyTinyBubbles cannot get value\n" );
		}
		else {

			////////////////////////////////////////////////////
			// read input values from handles
			////////////////////////////////////////////////////

			float input_val						= input_data.asFloat();
			MTime time_val						= time_data.asTime();
			MObject emitter_mesh_val			= emitter_mesh_data.asMesh();
			MString emitter_mesh_name_val		= emitter_mesh_name_data.asString();
			MString fluid_container_name_val	= fluid_container_name_data.asString();
			int emission_rate_val				= emission_rate_data.asInt();
			float scattering_frequency_val		= scattering_frequency_data.asFloat();
			float scattering_coefficient_val	= scattering_coefficient_data.asFloat();
			float breakup_frequency_val			= breakup_frequency_data.asFloat();
			float bubble_size_min_val			= bubble_size_min_data.asFloat();
			float bubble_size_max_val			= bubble_size_max_data.asFloat();


			// TODO: get fluid container attributes

			// TODO: get bubble emitter attributes



			// Get a handle to the output attribute.  This is similar to the
			// "inputValue" call above except that no dependency graph 
			// computation will be done as a result of this call.
			MDataHandle outputHandle = data.outputValue( ManyTinyBubbles::output );
			// This just copies the input value through to the output.  
			outputHandle.set( input_val );
			// Mark the destination plug as being clean.  This will prevent the
			// dependency graph from repeating this calculation until an input 
			// of this node changes.
			data.setClean( plug );

		}
	}
	else {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}



////////////////////////////////////////////////////
// creator(): exists to give Maya a way to create new objects of this type
//		returns a new object of this type
////////////////////////////////////////////////////
void* ManyTinyBubbles::creator()
{
	return new ManyTinyBubbles();
}


////////////////////////////////////////////////////
// initialize(): called to create and initialize all attributes and attribute dependencies for this node type;
//				 only called once when the node type is registered with Maya
////////////////////////////////////////////////////
MStatus ManyTinyBubbles::initialize()	
{
	////////////////////////////////////////////////////
	// create attributes
	////////////////////////////////////////////////////

	MFnNumericAttribute nAttr;
	MFnTypedAttribute	tAttr;
	MFnUnitAttribute	uAttr;
	MStatus				stat;


	////////////////////////////////////////////////////
	// create attributes
	////////////////////////////////////////////////////

	// setStorable: attribute will or will not be written to files when this type of node is stored
	// setKeyable: attribute is or is not keyable and will show up in the channel box

	ManyTinyBubbles::input = nAttr.create( "input", "in", MFnNumericData::kFloat, 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( false );

	ManyTinyBubbles::output = nAttr.create( "output", "out", MFnNumericData::kFloat, 0.0 );
	nAttr.setWritable( false );
	nAttr.setStorable( false );

	ManyTinyBubbles::time = uAttr.create( "time", "tm", MFnUnitAttribute::kTime, 0.0 );
	uAttr.setWritable( false );
 	uAttr.setStorable( false );

	ManyTinyBubbles::emitter_mesh = tAttr.create( "emitter_mesh", "em", MFnData::kMesh );
	//ManyTinyBubbles::emitter_mesh = tAttr.create( "emitter_mesh", "em", MFnMeshData::kMesh );
	tAttr.setStorable( true );
	tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::emitter_mesh_name = tAttr.create( "emitter_mesh_name", "emn", MFnData::kString );
 	tAttr.setStorable( true );
    tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::fluid_container_name = tAttr.create( "fluid_container_name", "fcn", MFnData::kString );
 	tAttr.setStorable( true );
    tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::emission_rate = nAttr.create( "emission_rate", "er", MFnNumericData::kInt, 1000 );
	nAttr.setSoftMax( 10000 );
	nAttr.setSoftMin( 1 );
	nAttr.setMin( 1 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::scattering_frequency = nAttr.create( "scattering_frequency", "sf", MFnNumericData::kFloat, 0.5 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( 0.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::scattering_coefficient = nAttr.create( "scattering_coefficient", "sc", MFnNumericData::kFloat, 0.0 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( -1.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( -1.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::breakup_frequency = nAttr.create( "breakup_frequency", "bf", MFnNumericData::kFloat, 0.5 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( 0.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::bubble_size_min = nAttr.create( "bubble_size_min", "smin", MFnNumericData::kFloat, 0.001 );
	nAttr.setSoftMax( 0.1 );
	nAttr.setSoftMin( 0.001 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::bubble_size_max = nAttr.create( "bubble_size_max", "smax", MFnNumericData::kFloat, 0.001 );
	nAttr.setSoftMax( 0.1 );
	nAttr.setSoftMin( 0.001 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );


	////////////////////////////////////////////////////
	// add attributes to node
	////////////////////////////////////////////////////

	stat = addAttribute( ManyTinyBubbles::input );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::time );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::emitter_mesh );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::emitter_mesh_name );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::fluid_container_name );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::emission_rate );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::scattering_frequency );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::scattering_coefficient );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::breakup_frequency );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::bubble_size_min );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::bubble_size_max );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }


	////////////////////////////////////////////////////
	// setup dependencies between attributes
	////////////////////////////////////////////////////

	stat = attributeAffects( ManyTinyBubbles::input, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::time, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::emitter_mesh, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::emitter_mesh_name, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::fluid_container_name, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::emission_rate, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::scattering_frequency, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::scattering_coefficient, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::breakup_frequency, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::bubble_size_min, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::bubble_size_max, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }

	return MS::kSuccess;
}