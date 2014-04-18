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

#include "Convenience.h"
#include "vec.h"


MTypeId ManyTinyBubbles::m_id( 0x70256 );


////////////////////////////////////////////////////
// node attributes
////////////////////////////////////////////////////

//MObject ManyTinyBubbles::input;
MObject ManyTinyBubbles::m_output;
MObject ManyTinyBubbles::m_time;
//MObject ManyTinyBubbles::m_emitter_mesh;
MObject ManyTinyBubbles::m_emitter_mesh_name;
MObject ManyTinyBubbles::m_fluid_container_name;
MObject ManyTinyBubbles::m_emission_rate;
MObject ManyTinyBubbles::m_scattering_frequency;
MObject ManyTinyBubbles::m_scattering_coefficient;
MObject ManyTinyBubbles::m_breakup_frequency;
MObject ManyTinyBubbles::m_bubble_size_min;
MObject ManyTinyBubbles::m_bubble_size_max;
MObject ManyTinyBubbles::m_step_size;


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
	// TODO: arrange output attribute in such a way that Maya recomputes it when I want it to

	// debug
	Convenience::printInScriptEditor( MString( "in compute()" ) );

	// debug
	//for ( unsigned int i = 1; i <= 10; ++i ) {
	//	testCode( i );
	//}

	MStatus returnStatus;


	////////////////////////////////////////////////////
	// check which output attribute we have been asked to compute
	////////////////////////////////////////////////////

	if ( plug == ManyTinyBubbles::m_output ) {

		// debug
		Convenience::printInScriptEditor( MString( "recomputing stuff" ) );


		////////////////////////////////////////////////////
		// get handles to input attribute we will need for computation
		////////////////////////////////////////////////////

		//MDataHandle input_data					= data.inputValue( ManyTinyBubbles::m_input, &returnStatus );
		MDataHandle time_data						= data.inputValue( ManyTinyBubbles::m_time, &returnStatus );
		//MDataHandle emitter_mesh_data				= data.inputValue( ManyTinyBubbles::m_emitter_mesh, &returnStatus );
		MDataHandle emitter_mesh_name_data			= data.inputValue( ManyTinyBubbles::m_emitter_mesh_name, &returnStatus );
		MDataHandle fluid_container_name_data		= data.inputValue( ManyTinyBubbles::m_fluid_container_name, &returnStatus );
		MDataHandle emission_rate_data				= data.inputValue( ManyTinyBubbles::m_emission_rate, &returnStatus );
		MDataHandle scattering_frequency_data		= data.inputValue( ManyTinyBubbles::m_scattering_frequency, &returnStatus );
		MDataHandle scattering_coefficient_data		= data.inputValue( ManyTinyBubbles::m_scattering_coefficient, &returnStatus );
		MDataHandle breakup_frequency_data			= data.inputValue( ManyTinyBubbles::m_breakup_frequency, &returnStatus );
		MDataHandle bubble_size_min_data			= data.inputValue( ManyTinyBubbles::m_bubble_size_min, &returnStatus );
		MDataHandle bubble_size_max_data			= data.inputValue( ManyTinyBubbles::m_bubble_size_max, &returnStatus );
		MDataHandle step_size_data					= data.inputValue( ManyTinyBubbles::m_step_size, &returnStatus );

		if ( returnStatus != MS::kSuccess ) {
			MGlobal::displayError( "Node ManyTinyBubbles cannot get value\n" );
		}
		else {

			////////////////////////////////////////////////////
			// read input values from handles
			////////////////////////////////////////////////////

			//float input_val					= input_data.asFloat();
			MTime time_val						= time_data.asTime();
			//MObject emitter_mesh_val			= emitter_mesh_data.asMesh();
			MString emitter_mesh_name_val		= emitter_mesh_name_data.asString();
			MString fluid_container_name_val	= fluid_container_name_data.asString();
			int emission_rate_val				= emission_rate_data.asInt();
			float scattering_frequency_val		= scattering_frequency_data.asFloat();
			float scattering_coefficient_val	= scattering_coefficient_data.asFloat();
			float breakup_frequency_val			= breakup_frequency_data.asFloat();
			float bubble_size_min_val			= bubble_size_min_data.asFloat();
			float bubble_size_max_val			= bubble_size_max_data.asFloat();
			float step_size_val					= step_size_data.asFloat();


			////////////////////////////////////////////////////
			// get fluid container attributes
			////////////////////////////////////////////////////

			MIntArray fluid_container_res_array = Convenience::getAttributeIntArray( fluid_container_name_val, MString( "resolution" ) );
			MDoubleArray fluid_container_dim_array = Convenience::getAttributeDoubleArray( fluid_container_name_val, MString( "dimensions" ) );

			// get fluid shape parent to retrieve translation attributes of fluid container
			MString fluid = Convenience::getParent( fluid_container_name_val );
			MDoubleArray fluid_container_translation_array = Convenience::getAttributeDoubleArray( fluid, MString( "translate" ) );

			m_fluid_container_res_x = fluid_container_res_array[VX];
			m_fluid_container_res_y = fluid_container_res_array[VY];
			m_fluid_container_res_z = fluid_container_res_array[VZ];

			m_fluid_container_dim_x = fluid_container_dim_array[VX];
			m_fluid_container_dim_y = fluid_container_dim_array[VY];
			m_fluid_container_dim_z = fluid_container_dim_array[VZ];

			m_fluid_container_trans_x = fluid_container_translation_array[VX];
			m_fluid_container_trans_y = fluid_container_translation_array[VY];
			m_fluid_container_trans_z = fluid_container_translation_array[VZ];

			m_fluid_container_cell_size_x = m_fluid_container_dim_x / m_fluid_container_res_x;
			m_fluid_container_cell_size_y = m_fluid_container_dim_y / m_fluid_container_res_y;
			m_fluid_container_cell_size_z = m_fluid_container_dim_z / m_fluid_container_res_z;


			////////////////////////////////////////////////////
			// get bubble emitter attributes
			////////////////////////////////////////////////////

			// TODO: implement this


			////////////////////////////////////////////////////
			// create bubbles
			////////////////////////////////////////////////////

			//createBubble(times, timeSteps, scatterFreqs, scatterCoefs, bubbleBreakFreqs, bubbleMinRadiuss, bubbleMaxRadiuss, newOutputData, plug, data );


			////////////////////////////////////////////////////
			// get a handle to the output attribute
			////////////////////////////////////////////////////

			MDataHandle output_data = data.outputValue( ManyTinyBubbles::m_output );
			
			// this just copies the input value through to the output
			//output_handle.set( input_val );
			output_data.set( ( int )time_val.as( MTime::kFilm ) );


			////////////////////////////////////////////////////
			// mark destination plug as clean to prevent dependency graph from repeating this calculation until an input of this node changes
			////////////////////////////////////////////////////

			data.setClean( plug );
		}
	}
	else {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}


////////////////////////////////////////////////////
// createBubbles()
////////////////////////////////////////////////////
MStatus ManyTinyBubbles::createBubbles( const MTime& time, float step_size, float scatter_freq,
										float scatter_coeff, float breakup_freq, float radius_min,
										float radius_max, MObject& out_data, const MPlug& plug,
										MDataBlock& block )
{
	// generate bubble radii based on radius_min and radius_max
	m_bubbles.setRadii( radius_min, radius_max );

	// ensure frame_num is not zero 
	int	frame_num = ( int )time.as( MTime::kFilm );
	if ( frame_num == 0 ) {
		frame_num = 1;
	}

	// delete all particles in Maya
	m_bubbles.deleteAllParticlesInMaya();


	// TODO: more


	return MS::kSuccess;
}


////////////////////////////////////////////////////
// creator(): exists to give Maya a way to create new objects of this type
//		returns a new object of this type
////////////////////////////////////////////////////
void* ManyTinyBubbles::creator()
{
	// debug
	Convenience::printInScriptEditor( MString( "in creator()" ) );

	return new ManyTinyBubbles();
}


////////////////////////////////////////////////////
// initialize(): called to create and initialize all attributes and attribute dependencies for this node type;
//				 only called once when the node type is registered with Maya
////////////////////////////////////////////////////
MStatus ManyTinyBubbles::initialize()	
{
	// debug
	Convenience::printInScriptEditor( MString( "in initialize()" ) );


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

	//ManyTinyBubbles::m_input = nAttr.create( "input", "in", MFnNumericData::kFloat, 0.0 );
	//nAttr.setStorable( true );
	//nAttr.setKeyable( false );

	ManyTinyBubbles::m_output = nAttr.create( "output", "out", MFnNumericData::kFloat, 0.0 );
	nAttr.setWritable( false );
	nAttr.setStorable( false );
	nAttr.setKeyable( false );

	ManyTinyBubbles::m_time = uAttr.create( "time", "tm", MFnUnitAttribute::kTime, 0.0 );
	uAttr.setWritable( false );
 	uAttr.setStorable( false );

	//ManyTinyBubbles::m_emitter_mesh = tAttr.create( "emitter_mesh", "em", MFnData::kMesh );
	//ManyTinyBubbles::m_emitter_mesh = tAttr.create( "emitter_mesh", "em", MFnMeshData::kMesh );
	//tAttr.setStorable( true );
	//tAttr.setWritable( false );
	//tAttr.setStorable( false );

	ManyTinyBubbles::m_emitter_mesh_name = tAttr.create( "emitter_mesh_name", "emn", MFnData::kString );
 	tAttr.setStorable( true );
    tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::m_fluid_container_name = tAttr.create( "fluid_container_name", "fcn", MFnData::kString );
 	tAttr.setStorable( true );
    tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::m_emission_rate = nAttr.create( "emission_rate", "er", MFnNumericData::kInt, 1000 );
	nAttr.setSoftMax( 10000 );
	nAttr.setSoftMin( 1 );
	nAttr.setMin( 1 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::m_scattering_frequency = nAttr.create( "scattering_frequency", "sf", MFnNumericData::kFloat, 0.5 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( 0.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::m_scattering_coefficient = nAttr.create( "scattering_coefficient", "sc", MFnNumericData::kFloat, 0.0 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( -1.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( -1.0 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::m_breakup_frequency = nAttr.create( "breakup_frequency", "bf", MFnNumericData::kFloat, 0.5 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( 0.0 );
	nAttr.setMax( 1.0 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::m_bubble_size_min = nAttr.create( "bubble_size_min", "smin", MFnNumericData::kFloat, 0.001 );
	nAttr.setSoftMax( 0.1 );
	nAttr.setSoftMin( 0.001 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::m_bubble_size_max = nAttr.create( "bubble_size_max", "smax", MFnNumericData::kFloat, 0.001 );
	nAttr.setSoftMax( 0.1 );
	nAttr.setSoftMin( 0.001 );
	nAttr.setMin( 0.0 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::m_step_size = nAttr.create( "simulation_step_size", "sss", MFnNumericData::kFloat, 0.1 );
	nAttr.setSoftMax( 1.0 );
	nAttr.setSoftMin( 0.001 );
	nAttr.setMin( 0.001 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );


	////////////////////////////////////////////////////
	// add attributes to node
	////////////////////////////////////////////////////

	//stat = addAttribute( ManyTinyBubbles::m_input );
	//	if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_time );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	//stat = addAttribute( ManyTinyBubbles::m_emitter_mesh );
	//	if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_emitter_mesh_name );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_fluid_container_name );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_emission_rate );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_scattering_frequency );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_scattering_coefficient );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_breakup_frequency );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_bubble_size_min );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_bubble_size_max );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_step_size );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }


	////////////////////////////////////////////////////
	// setup dependencies between attributes
	////////////////////////////////////////////////////

	//stat = attributeAffects( ManyTinyBubbles::m_input, ManyTinyBubbles::m_output );
	//	if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_time, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	//stat = attributeAffects( ManyTinyBubbles::m_emitter_mesh, ManyTinyBubbles::m_output );
	//	if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_emitter_mesh_name, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_fluid_container_name, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_emission_rate, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_scattering_frequency, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_scattering_coefficient, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_breakup_frequency, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_bubble_size_min, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_bubble_size_max, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_step_size, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }

	return MS::kSuccess;
}


////////////////////////////////////////////////////
// debug
////////////////////////////////////////////////////
void ManyTinyBubbles::testCode( const unsigned int& num ) const
{
	std::string cmd = "particleExists bubbleParticle";
	Convenience::appendNumToStdString( cmd, num );
	cmd += ";";

	// debug
	Convenience::printInScriptEditor( cmd );

	//int particle_exists;
	//	MGlobal::executeCommand( Convenience::convertStdStringToMString( cmd ), particle_exists );
}