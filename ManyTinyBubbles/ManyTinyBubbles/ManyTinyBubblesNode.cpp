//
// Copyright 2014 (C) Danny Rerucha and Wei-Chien Tu
// 
// File: ManyTinyBubblesNode.cpp
//
// Dependency Graph Node: ManyTinyBubbles
//
// Authors: Danny Rerucha and Wei-Chien Tu
//


// TODO: add if defines to all classes
// TODO: change all floats to doubles b/c I think they're faster on x86 architecture


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


////////////////////////////////////////////////////
// constants
////////////////////////////////////////////////////

const double GAS_DENSITY = 100.0;


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

ManyTinyBubbles::ManyTinyBubbles()
{
	m_current_frame = 0;
}

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

			// get fluid container attributes and store in m_fluid_container
			m_fluid_container.init( fluid_container_name_val );

			// store emitter in m_bubbles
			m_bubbles.init( scattering_frequency_val,
							scattering_coefficient_val,
							breakup_frequency_val,
							bubble_size_min_val,
							bubble_size_max_val );


			// TODO: get bubble emitter attributes


			////////////////////////////////////////////////////
			// create bubbles
			////////////////////////////////////////////////////

			createBubbles( time_val,
						   step_size_val );


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
MStatus ManyTinyBubbles::createBubbles( const MTime& time,
										const float& step_size )
{
	// ensure frame_num is not zero 
	unsigned int frame_num = ( unsigned int )time.as( MTime::kFilm );
	if ( frame_num == 0 ) {
		frame_num = 1;
	}

	// if the user wants a previously simulated frame, then we must run the simulation from the very beginning
	if ( frame_num < m_current_frame ) {
		reset();
		m_current_frame = 0;
	}

	// delete all particles in Maya
	m_bubbles.deleteAllParticlesInMaya();

	// simulate ( frame_num - m_current_frame ) frames starting from current m_current_frame
	for ( unsigned int i = 0; i < frame_num - m_current_frame; ++i ) {
		computeFractionField();

		// TODO: set fluid density
		//setDensity( "fluid1" );

		advectParticles( step_size );

		//computeBubbleGenerationPosFromMesh( LEVEL_SET_RES, 100 );
		//generateBubbles();
	}

	m_current_frame = frame_num;






	// TODO: create particle group in Maya

	// TODO: set bubble size in Maya

	// TODO: maybe create a dummy mesh if the program doesn't work without one?


	return MS::kSuccess;
}


////////////////////////////////////////////////////
// computeFractionField()
////////////////////////////////////////////////////
void ManyTinyBubbles::computeFractionField()
{
	// set fraction field value at each voxel to 1.0f
	m_fluid_container.resetFractionField();

	// get positions for every bubble present in the simulation
	std::vector<std::vector<vec3>> bubble_pos_list = m_bubbles.getPosList();

	// iterate through the bubble radius groups
	unsigned int radius_index = 0;
	for ( std::vector<std::vector<vec3>>::iterator outer_it = bubble_pos_list.begin() ; outer_it != bubble_pos_list.end(); ++outer_it ) {
		std::vector<vec3> bubble_pos_sublist = *outer_it;

		// iterate through the list of vec3s in the current radius group
		for ( std::vector<vec3>::iterator inner_it = bubble_pos_sublist.begin() ; inner_it != bubble_pos_sublist.end(); ++inner_it ) {
			vec3 bubble_pos = *inner_it;

			// reduce fraction field of the voxel the current bubble is inside
			m_fluid_container.reduceFractionFieldOfVoxelAtPos( bubble_pos, m_bubbles.getRadiusAtIndex( radius_index ) );
		}

		++radius_index;
	}
}


////////////////////////////////////////////////////
// simulate bubbles
////////////////////////////////////////////////////
void ManyTinyBubbles::advectParticles( const float& dt )
{
	// update velocity field to match the Maya fluid
	m_fluid_container.updateVelocityField();

	// get positions for every bubble present in the simulation
	std::vector<std::vector<vec3>> bubble_pos_list = m_bubbles.getPosList();

	// iterate through the bubble radius groups
	for ( std::vector<std::vector<vec3>>::iterator outer_it = bubble_pos_list.begin() ; outer_it != bubble_pos_list.end(); ++outer_it ) {
		std::vector<vec3> bubble_pos_sublist = *outer_it;

		// iterate through the list of vec3s in the current radius group
		for ( std::vector<vec3>::iterator inner_it = bubble_pos_sublist.begin() ; inner_it != bubble_pos_sublist.end(); ++inner_it ) {
			vec3 bubble_pos = *inner_it;

			// TODO: get actual bubble velocity here instead of just the velocity of the voxel the bubble is in
			vec3 bubble_vel = m_fluid_container.getVelocityOfVoxelAtPos( bubble_pos );

			double scattering_probability = computeScatteringProbabilityOfBubble( bubble_vel, bubble_pos );

			// bubble will scatter, so alter bubble direction
			if ( scattering_probability > Convenience::generateRandomDoubleBetweenZeroAndOneInclusive() ) {
				double altered_angle = computeAlteredAngle();

				// altered angle is not zero, so update bubble velocity direction
 				if ( altered_angle < -DBL_EPSILON || altered_angle > DBL_EPSILON ) {

					// update bubble velocity due to scattering
					bubble_vel = updateBubbleVelocity( bubble_vel, altered_angle );
				}
			}

			// TODO: update bubble position using the altered velocity
			// TODO: update all bubble positions even if their velocity was not altered, probably

			// TODO: check if new bubble position escapes fluid container


		}
	}
}


////////////////////////////////////////////////////
// compute scattering probability function, s(x) in paper, [0, 1]
////////////////////////////////////////////////////
double ManyTinyBubbles::computeScatteringProbabilityOfBubble( const vec3& bubble_vel,
															  const vec3& bubble_pos ) const
{
	double voxel_vel_magnitude = bubble_vel.Length();
	double fraction_field = m_fluid_container.getFractionFieldOfVoxelAtPos( bubble_pos );

	return m_bubbles.getScatteringFrequency() * GAS_DENSITY * ( 1.0 - fraction_field ) * voxel_vel_magnitude * voxel_vel_magnitude;
}


////////////////////////////////////////////////////
// compute bubble's new direction
////////////////////////////////////////////////////
double ManyTinyBubbles::computeAlteredAngle() const
{
	double uniform_random_num = Convenience::generateRandomDoubleBetweenZeroAndOneInclusive();
	double numerator = 2.0 * uniform_random_num + m_bubbles.getScatteringCoefficient() - 1.0;
	double denominator = 2.0 * m_bubbles.getScatteringCoefficient() * uniform_random_num - m_bubbles.getScatteringCoefficient() + 1.0;

	// compute theta, altered angle
	double theta;
	if ( denominator > -DBL_EPSILON && denominator < DBL_EPSILON  ) {
		// if denominator == zero
		theta = M_PI / 2.0;
	}
	else {
		theta = acos( numerator / denominator );
	}

	return theta;
}


////////////////////////////////////////////////////
// updateBubbleVelocity()
////////////////////////////////////////////////////
vec3 ManyTinyBubbles::updateBubbleVelocity( const vec3&		old_vel,
											const double&	altered_dir ) const
{
	vec3 new_vel;

	// create a random axis to rotate around using quaternions
	// random axis must be perpendicular to velocity direction
	// so, the dot product of the random axis and the velocity direction should equal zero
	// or, rotation_axis_x * new_vel[VX] + rotation_axis_y * new_vel[VY] + rotation_axis_z * new_vel[VZ] == 0
	// we randomly provide the first two values and use the above equation to compute the third
	// we make careful considerations to avoid divide by zero

	// TODO: ask why these variables are initialized this way

	vec3 rotation_axis( 1.0, 0.0, 0.0 );
					
	// if velocity[VZ] == 0, then the computation of rotation_axis_z will divide by zero
	if ( old_vel[VZ] > -DBL_EPSILON && old_vel[VZ] < DBL_EPSILON ) {

		// generate random numbers for x and y components [-1, 1] so that they are not both 0
		do {
			rotation_axis[VX] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
			rotation_axis[VY] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
		} while ( ( rotation_axis[VX] > -DBL_EPSILON && rotation_axis[VX] < DBL_EPSILON ) &&
				  ( rotation_axis[VY] > -DBL_EPSILON && rotation_axis[VY] < DBL_EPSILON ) );

		rotation_axis[VZ] = -1.0 * ( rotation_axis[VX] * old_vel[VX] + rotation_axis[VY] * old_vel[VY] ) / old_vel[VZ];
	}
	else if ( old_vel[VY] > -DBL_EPSILON && old_vel[VY] < DBL_EPSILON ) {

		// generate random numbers for x and z components [-1, 1] so that they are not both 0
		do {
			rotation_axis[VX] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
			rotation_axis[VZ] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
		} while ( ( rotation_axis[VX] > -DBL_EPSILON && rotation_axis[VX] < DBL_EPSILON ) &&
				  ( rotation_axis[VZ] > -DBL_EPSILON && rotation_axis[VZ] < DBL_EPSILON ) );

		rotation_axis[VZ] = -1.0 * ( rotation_axis[VX] * old_vel[VX] + rotation_axis[VZ] * old_vel[VZ] ) / old_vel[VY];
	}
	else if ( old_vel[VX] > -DBL_EPSILON && old_vel[VX] < DBL_EPSILON ) {

		// generate random numbers for y and z components [-1, 1] so that they are not both 0
		do {
			rotation_axis[VY] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
			rotation_axis[VZ] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
		} while ( ( rotation_axis[VY] > -DBL_EPSILON && rotation_axis[VY] < DBL_EPSILON ) &&
				  ( rotation_axis[VZ] > -DBL_EPSILON && rotation_axis[VZ] < DBL_EPSILON ) );

		rotation_axis[VX] = -1.0 * ( rotation_axis[VY] * old_vel[VY] + rotation_axis[VZ] * old_vel[VZ] ) / old_vel[VX];

	}

	// TODO: check for potential divide by zero before calling Normalize()

	// make unit length
	rotation_axis.Normalize();

	// use quaternion rotation to rotate velocity
	double half_altered_dir = altered_dir / 2.0;
	double sin_half_altered_dir = sin( half_altered_dir );
	double ss = cos( half_altered_dir ); 
	double xx = sin_half_altered_dir * rotation_axis[VX];
	double yy = sin_half_altered_dir * rotation_axis[VY];
	double zz = sin_half_altered_dir * rotation_axis[VZ];

	// convert quaternion to rotation matrix and multiply by original velocity vector
	new_vel[VX] = ( 1.0 - 2.0 * yy * yy - 2.0 * zz * zz ) * old_vel[VX] + ( 2.0 * xx * yy - 2.0 * ss * zz ) * old_vel[VY] + (2.0 * xx * zz + 2.0 * ss * yy ) * old_vel[VZ];
	new_vel[VY] = ( 2.0 * xx * yy + 2.0 * ss * zz ) * old_vel[VX] + ( 1.0 - 2.0 * xx * xx - 2.0 * zz * zz ) * old_vel[VY] + ( 2.0 * yy * zz - 2.0 * ss * xx ) * old_vel[VZ];
	new_vel[VZ] = ( 2.0 * xx * zz - 2.0 * ss * yy ) * old_vel[VX] + ( 2.0 * yy * zz + 2.0 * ss * xx ) * old_vel[VY] + ( 1.0 - 2.0 * xx * xx - 2.0 * yy * yy ) * old_vel[VZ];

	return new_vel;
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


void ManyTinyBubbles::reset()
{
	m_bubbles.reset();

	// TODO: set fraction field to some default value
	// TODO: allow fraction field to handle voxels that aren't squares
	//mFractionField.initialize(1.0f, m_fluid_container_cell_size_x, m_fluid_container_dim_x, m_fluid_container_dim_y, m_fluid_container_dim_z); // set default fraction field 
}


////////////////////////////////////////////////////
// debug
////////////////////////////////////////////////////
void ManyTinyBubbles::testCode( const MString& str  ) const
{
	//MDoubleArray velocity_field = Convenience::getVelocityFieldFromMayaFluid( str );

	//std::string to_print = "velocity_field size: ";
	//Convenience::appendNumToStdString( to_print, velocity_field.length() );
	//Convenience::printInScriptEditor( to_print );
}