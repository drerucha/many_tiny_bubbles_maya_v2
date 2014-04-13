//
// Copyright (C) CIS660
// 
// File: ManyTinyBubblesNode.cpp
//
// Dependency Graph Node: ManyTinyBubbles
//
// Author: Maya Plug-in Wizard 2.0
//

#include "ManyTinyBubblesNode.h"

#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

#include <maya/MGlobal.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnTypedAttribute.h>


MTypeId     ManyTinyBubbles::id( 0x70256 );

// node attributes

MObject ManyTinyBubbles::input;		// test
MObject ManyTinyBubbles::output;	// test

MObject ManyTinyBubbles::time;
MObject ManyTinyBubbles::input_mesh;
MObject ManyTinyBubbles::emission_rate;
MObject ManyTinyBubbles::scattering_frequency;
MObject ManyTinyBubbles::scattering_coefficient;
MObject ManyTinyBubbles::breakup_frequency;
MObject ManyTinyBubbles::bubble_size_min;
MObject ManyTinyBubbles::bubble_size_max;


// constructor/destructor
ManyTinyBubbles::ManyTinyBubbles() {}
ManyTinyBubbles::~ManyTinyBubbles() {}


//
//	Description:
//		This method computes the value of the given output plug based
//		on the values of the input attributes.
//
//	Arguments:
//		plug - the plug to compute
//		data - object that provides access to the attributes for this node
//
MStatus ManyTinyBubbles::compute( const MPlug& plug, MDataBlock& data )
{
	MStatus returnStatus;
 
	// Check which output attribute we have been asked to compute.  If this 
	// node doesn't know how to compute it, we must return 
	// MS::kUnknownParameter.
	if ( plug == output ) {
		// Get a handle to the input attribute that we will need for the
		// computation.  If the value is being supplied via a connection 
		// in the dependency graph, then this call will cause all upstream  
		// connections to be evaluated so that the correct value is supplied.
		MDataHandle inputData = data.inputValue( input, &returnStatus );

		if ( returnStatus != MS::kSuccess ) {
			MGlobal::displayError( "Node ManyTinyBubbles cannot get value\n" );
		}
		else {
			// Read the input value from the handle.
			float result = inputData.asFloat();

			// Get a handle to the output attribute.  This is similar to the
			// "inputValue" call above except that no dependency graph 
			// computation will be done as a result of this call.
			MDataHandle outputHandle = data.outputValue( ManyTinyBubbles::output );
			// This just copies the input value through to the output.  
			outputHandle.set( result );
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


//
//	Description:
//		this method exists to give Maya a way to create new objects
//      of this type. 
//
//	Return Value:
//		a new object of this type
//
void* ManyTinyBubbles::creator()
{
	return new ManyTinyBubbles();
}


//
//	Description:
//		This method is called to create and initialize all of the attributes
//      and attribute dependencies for this node type.  This is only called 
//		once when the node type is registered with Maya.
//
//	Return Values:
//		MS::kSuccess
//		MS::kFailure
//	
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
	// test attributes
	////////////////////////////////////////////////////

	ManyTinyBubbles::input = nAttr.create( "input", "in", MFnNumericData::kFloat, 0.0 );
	// Attribute will be written to files when this type of node is stored
 	nAttr.setStorable( true );
	// Attribute is keyable and will show up in the channel box
 	nAttr.setKeyable( false );

	ManyTinyBubbles::output = nAttr.create( "output", "out", MFnNumericData::kFloat, 0.0 );
	// Attribute is read-only because it is an output attribute
	nAttr.setWritable( false );
	// Attribute will not be written to files when this type of node is stored
	nAttr.setStorable( false );


	////////////////////////////////////////////////////
	// my attributes
	////////////////////////////////////////////////////

	ManyTinyBubbles::time = uAttr.create( "time", "tm", MFnUnitAttribute::kTime, 0.0 );
	uAttr.setWritable( false );
 	uAttr.setStorable( false );

	//ManyTinyBubbles::input_mesh = tAttr.create( "mesh", "m", MFnData::kMesh );
	ManyTinyBubbles::input_mesh = tAttr.create( "mesh", "m", MFnMeshData::kMesh );
 	tAttr.setStorable( true );
    tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::emission_rate = nAttr.create( "rate", "r", MFnNumericData::kInt, 1000 );
	nAttr.setSoftMax( 10000 );
	nAttr.setSoftMin( 1 );
	nAttr.setMin( 1 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::scattering_frequency = nAttr.create( "scattfreq", "sf", MFnNumericData::kFloat, 0.5 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( 0.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::scattering_coefficient = nAttr.create( "scattcoeff", "sc", MFnNumericData::kFloat, 0.0 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( -1.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( -1.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::breakup_frequency = nAttr.create( "breakup", "brk", MFnNumericData::kFloat, 0.5 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( 0.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::bubble_size_min = nAttr.create( "sizemin", "min", MFnNumericData::kFloat, 0.001 );
	nAttr.setSoftMax( 0.1 );
	nAttr.setSoftMin( 0.001 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::bubble_size_max = nAttr.create( "sizemax", "max", MFnNumericData::kFloat, 0.001 );
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
	stat = addAttribute( ManyTinyBubbles::input_mesh );
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

	// Set up a dependency between the input and the output.  This will cause
	// the output to be marked dirty when the input changes.  The output will
	// then be recomputed the next time the value of the output is requested.
	stat = attributeAffects( ManyTinyBubbles::input, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }

	stat = attributeAffects( ManyTinyBubbles::time, ManyTinyBubbles::output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::input_mesh, ManyTinyBubbles::output );
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


	////////////////////////////////////////////////////
	// return value
	////////////////////////////////////////////////////

	return MS::kSuccess;
}