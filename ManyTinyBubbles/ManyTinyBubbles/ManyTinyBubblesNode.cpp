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
// TODO: alter node to allow multiple bubble emitter sources
// TODO: index vectors using [] instead of at() b/c it seems faster


#include "ManyTinyBubblesNode.h"

#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MGlobal.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MTime.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>

#include "Convenience.h"
#include "vec.h"
#include "GlobalState.h"


////////////////////////////////////////////////////
// constants
////////////////////////////////////////////////////

const double GAS_DENSITY = 100.0;
const std::string MAYA_PARTICLE_NAME = "bubbleParticle";


////////////////////////////////////////////////////
// node attributes
////////////////////////////////////////////////////

MTypeId ManyTinyBubbles::m_id( 0x70256 );

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

MObject ManyTinyBubbles::m_fluid_transform_name;
MObject ManyTinyBubbles::m_node_name;
MObject ManyTinyBubbles::m_melting_rate;
MObject ManyTinyBubbles::m_level_set_resolution;


////////////////////////////////////////////////////
// constructor/destructor
////////////////////////////////////////////////////

ManyTinyBubbles::ManyTinyBubbles()
{
	m_current_frame = 0;
}

ManyTinyBubbles::~ManyTinyBubbles()
{

}


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

	if ( plug == ManyTinyBubbles::m_output ) {

		////////////////////////////////////////////////////
		// get handles to input attribute we will need for computation
		////////////////////////////////////////////////////

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

		MDataHandle fluid_transform_name_data		= data.inputValue( ManyTinyBubbles::m_fluid_transform_name, &returnStatus );
		MDataHandle node_name_data					= data.inputValue( ManyTinyBubbles::m_node_name, &returnStatus );
		MDataHandle melting_rate_data				= data.inputValue( ManyTinyBubbles::m_melting_rate, &returnStatus );
		MDataHandle level_set_resolution_data		= data.inputValue( ManyTinyBubbles::m_level_set_resolution, &returnStatus );

		if ( returnStatus != MS::kSuccess ) {
			MGlobal::displayError( "Node ManyTinyBubbles cannot get value\n" );
		}
		else {

			////////////////////////////////////////////////////
			// read input values from handles
			////////////////////////////////////////////////////

			MTime time_val							= time_data.asTime();
			//MObject emitter_mesh_val				= emitter_mesh_data.asMesh();
			MString emitter_mesh_name_val			= emitter_mesh_name_data.asString();
			MString fluid_container_name_val		= fluid_container_name_data.asString();
			int emission_rate_val					= emission_rate_data.asInt();
			float scattering_frequency_val			= scattering_frequency_data.asFloat();
			float scattering_coefficient_val		= scattering_coefficient_data.asFloat();
			float breakup_frequency_val				= breakup_frequency_data.asFloat();
			float bubble_size_min_val				= bubble_size_min_data.asFloat();
			float bubble_size_max_val				= bubble_size_max_data.asFloat();
			float step_size_val						= step_size_data.asFloat();

			MString fluid_transform_name_val		= fluid_transform_name_data.asString();
			MString node_name_val					= node_name_data.asString();
			unsigned int emitter_level_set_res_val	= level_set_resolution_data.asInt();
			unsigned int emitter_melting_rate_val	= melting_rate_data.asInt();

			
			////////////////////////////////////////////////////
			// store fluid tranform name
			////////////////////////////////////////////////////

			// trim first char b/c we want something like "fluid1", not "|fluid1"
			std::string fluid_transform_name = Convenience::convertMStringToStdString( fluid_transform_name_val );
			fluid_transform_name = fluid_transform_name.substr( 1 );
			GlobalState::storeFluidTransformName( fluid_transform_name );


			////////////////////////////////////////////////////
			// initialize data structures
			////////////////////////////////////////////////////

			// TODO: only re-init these data structures if their respective data has changed

			// get fluid container attributes and store in m_fluid_container
			m_fluid_container.init( fluid_container_name_val );

			// store attributes governing bubble behavior in m_bubbles
			m_bubbles.init( scattering_frequency_val,
							scattering_coefficient_val,
							breakup_frequency_val,
							bubble_size_min_val,
							bubble_size_max_val );

			// store emitter name and infor (vertex list, face list) in GlobalState class
			// check if info for emitter is already stored to reduce rendundant computations
			if ( emitter_mesh_name_val.length() != 0 &&
				 GlobalState::objectExists( Convenience::convertMStringToStdString( emitter_mesh_name_val ) ) == false )
			{
				GlobalState::setSelectedObject( Convenience::convertMStringToStdString( emitter_mesh_name_val ) );

				int exists;
				MGlobal::executeCommand( "objExists " + emitter_mesh_name_val, exists );

				if ( exists ) {
					storeMeshInfoByName( emitter_mesh_name_val );
				}
			}

			m_emitter.init( emission_rate_val, bubble_size_min_val,
							emitter_level_set_res_val, emitter_melting_rate_val );


			////////////////////////////////////////////////////
			// get a handle to the output attribute
			////////////////////////////////////////////////////

			MDataHandle output_data = data.outputValue( ManyTinyBubbles::m_output );


			////////////////////////////////////////////////////
			// perform bubble simulation
			////////////////////////////////////////////////////

			simulationLoop( time_val, step_size_val );
			//simulationLoop( time_val, step_size_val, new_output_data,
			//				plug, data );


			////////////////////////////////////////////////////
			// create dummy cube mesh for output
			////////////////////////////////////////////////////

			MFnMeshData data_creator;
			MObject new_output_data = data_creator.create( &returnStatus );

			MFnMesh	meshFS;
			float cubeSize = 0.000001f ;
			MFloatPointArray points;
			const int numFaces = 6;
			int numVertices	= 8;
			const int numFaceConnects = 24;

			MFloatPoint vtx_1( -cubeSize, -cubeSize, -cubeSize );
			MFloatPoint vtx_2(  cubeSize, -cubeSize, -cubeSize );
			MFloatPoint vtx_3(  cubeSize, -cubeSize,  cubeSize );
			MFloatPoint vtx_4( -cubeSize, -cubeSize,  cubeSize );
			MFloatPoint vtx_5( -cubeSize,  cubeSize, -cubeSize );
			MFloatPoint vtx_6( -cubeSize,  cubeSize,  cubeSize );
			MFloatPoint vtx_7(  cubeSize,  cubeSize,  cubeSize );
			MFloatPoint vtx_8(  cubeSize,  cubeSize, -cubeSize );
			points.append( vtx_1 );
			points.append( vtx_2 );
			points.append( vtx_3 );
			points.append( vtx_4 );
			points.append( vtx_5 );
			points.append( vtx_6 );
			points.append( vtx_7 );
			points.append( vtx_8 );
			int face_counts[numFaces] = { 4, 4, 4, 4, 4, 4 };
			MIntArray faceCounts( face_counts, numFaces );
			int face_connects[ numFaceConnects ] = { 0, 1, 2, 3, 4, 5, 6, 7, 3, 2, 6, 5, 0, 3, 5, 4, 0, 4, 7, 1, 1, 7, 6, 2 };
			MIntArray faceConnects( face_connects, numFaceConnects );
			MObject newMesh = meshFS.create( numVertices, numFaces, points, faceCounts, faceConnects, new_output_data );

			//output_data.set( ( int )time_val.as( MTime::kFilm ) );
			output_data.set( new_output_data );


			////////////////////////////////////////////////////
			// mark destination plug as clean to prevent dependency graph from repeating this calculation until an input of this node changes
			////////////////////////////////////////////////////

			data.setClean( plug );

			// select Many Tiny Bubbles node
			MGlobal::executeCommand( "select -replace " + node_name_val );
		}
	}
	else {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}


////////////////////////////////////////////////////
// simulationLoop()
////////////////////////////////////////////////////
MStatus	ManyTinyBubbles::simulationLoop( const MTime& time, const float& step_size )
{
	// cleanup before simulation
	unsigned int requested_frame = simulationSetup( time );

	// loop from current frame to requested frame
	for ( unsigned int i = 0; i < requested_frame - m_current_frame; ++i ) {

		// TODO: do not use hardcoded name here
		// TODO: move hardcoded 1000 value below to the top of this file for easier access

		// if polySurface1 exists, then Maya fluid has been converted to polygon mesh and we should use level set method to determine fluid surface boundary
		int exists;
		MGlobal::executeCommand( "objExists polySurface1", exists );
		if ( exists == 1 ) {

			// compute and store vertex list and face list of fluid for level set method
			m_fluid_container.prepareFluidMeshForLevelSetMethod( "polySurface1" );

			// for many particles, compute signed distance function of entire fluid which is used later to determine which particles lie near the fluid surface
			if ( m_bubbles.getTotalBubbleNumber() > 1000 ) {
				m_fluid_container.updateLevelSetSignedDistanceFunction( "polySurface1" );
			}
		}

		// update bubble positions using explicit Euler integration
		m_bubbles.updateBubblePositions( step_size );

		// update density field to match the Maya fluid
		m_fluid_container.updateDensityField();

		// remove bubbles that "escape" fluid container
		deleteEscapedBubbles();

		// generate new bubbles on emitter surface
		generateMoreBubblesFromEmitter();

		// update bubble velocities taking their current velocities, the fluid velocity field, and their scattering properties into account
		updateBubbleVelocities();

		// bubbles randomly break apart into two smaller bubbles
		breakupBubbles();
	}

	// render bubbles as particles in Maya
	m_bubbles.createMayaParticlesWithName( MAYA_PARTICLE_NAME );		// create particle groups
	m_bubbles.setRadiiForMayaParticlesWithName( MAYA_PARTICLE_NAME );	// set particle radii

	// update m_current_frame
	m_current_frame = requested_frame;

	return MS::kSuccess;
}


////////////////////////////////////////////////////
// simulationSetup()
////////////////////////////////////////////////////
unsigned int ManyTinyBubbles::simulationSetup( const MTime& time )
{
	// ensure we simulate at least one frame to visualize results
	unsigned int requested_frame = ( unsigned int )time.as( MTime::kFilm );
	if ( requested_frame == 0 ) {
		requested_frame = 1;
	}

	// if the user wants to visualize a previously simulated frame, then we must re-run the simulation from the beginning
	if ( requested_frame < m_current_frame ) {
		reset();
		m_current_frame = 0;

		// remove current mesh from scene
		m_emitter.deleteEmitterMeshesFromScene();

		// create obj file and import it from stored face list data and vertex list data
		m_emitter.createObjFileFromStoredMeshData();
	}

	// delete all rendered particles in Maya b/c we recreate them each simulation
	m_bubbles.deleteAllParticlesInMaya();

	return requested_frame;
}


////////////////////////////////////////////////////
// reset()
////////////////////////////////////////////////////
void ManyTinyBubbles::reset()
{
	m_bubbles.reset();

	// set fraction field value at each voxel to 1.0f
	m_fluid_container.resetFractionField();
}


////////////////////////////////////////////////////
// deleteEscapedBubbles()
////////////////////////////////////////////////////
void ManyTinyBubbles::deleteEscapedBubbles()
{
	// TODO: remove hardcoded name "polySurface1"
	// TODO: move hardcoded value "1000" to the top of this file

	std::vector<BubbleLocator> bubbles_to_remove;

	// loop through radius groups
	unsigned int num_radius_groups = m_bubbles.getNumRadiusGroups();
	for ( unsigned int radius_group_index = 0; radius_group_index < num_radius_groups; ++radius_group_index ) {

		// loop through bubbles in the current radius group
		unsigned int num_bubbles_in_list = m_bubbles.getNumBubblesInListWithIndex( radius_group_index );
		for ( unsigned int bubble_index = 0; bubble_index < num_bubbles_in_list; ++bubble_index ) {

			// get bubble position
			vec3 bubble_pos = m_bubbles.getPosAtIndex( radius_group_index, bubble_index );

			// mark the bubble for deletion if it has "escaped" the fluid container
			//if ( m_fluid_container.posIsOutsideFluidContainer( bubble_pos ) ) {
			//	bubbles_to_remove.push_back( BubbleLocator( radius_group_index, bubble_index ) );
			//}

			int exists;
			MGlobal::executeCommand( "objExists polySurface1", exists );

			// if polySurface1 exists, then fluid has been converted into a polygon mesh so we use the level set method to determine fluid surface boundary
			if ( exists == 1 ) {

				// if too many bubbles, then filter bubbles so we only concern ourselves with bubbles near the fluid's surface
				if ( m_bubbles.getTotalBubbleNumber() > 1000 ) {

					// bubbles are outside fluid surface boundary
					if ( m_fluid_container.posIsUnderBoundary( bubble_pos ) == 1 ) {
						bubbles_to_remove.push_back( BubbleLocator( radius_group_index, bubble_index ) );
					}
					// bubbles are located near the fluid surface boundary
					else if ( m_fluid_container.posIsUnderBoundary( bubble_pos ) == 0 ) {

						// use finer level set method to determine if bubbles should be removed
						if ( m_fluid_container.isPosOutsideFluidViaLevelSet( bubble_pos ) ) {
							bubbles_to_remove.push_back( BubbleLocator( radius_group_index, bubble_index ) );
						}
					}
				}
				// if there aren't too many bubbles, directly use the finer level set method to filter bubbles
				else {
					if ( m_fluid_container.isPosOutsideFluidViaLevelSet( bubble_pos ) ) {
						bubbles_to_remove.push_back( BubbleLocator( radius_group_index, bubble_index ) );
					}
				}
			}
			else
			{
				if ( m_fluid_container.isPosOutsideFluidViaDensity( bubble_pos ) ) {
					bubbles_to_remove.push_back( BubbleLocator( radius_group_index, bubble_index ) );
				}
			}
		}
	}

	deleteBubblesInList( bubbles_to_remove );
}


////////////////////////////////////////////////////
// deleteBubblesInList()
////////////////////////////////////////////////////
void ManyTinyBubbles::deleteBubblesInList( std::vector<BubbleLocator> bubbles_to_remove )
{
	// TODO: fix this hack

	unsigned int counter = 0;
	unsigned int radius_group_index = 0;
	bool firstTimeFlag = true;

	// delete bubbles that were marked for deletion
	for ( std::vector<BubbleLocator>::iterator it = bubbles_to_remove.begin() ; it != bubbles_to_remove.end(); ++it ) {

		if ( it->radius_group != radius_group_index ) {
			radius_group_index = it->radius_group;
			counter = 0;
		}

		m_bubbles.removeBubbleAtIndex( it->radius_group, it->list_index - counter );

		++counter;
	}
}


////////////////////////////////////////////////////
// generateMoreBubblesFromEmitter()
////////////////////////////////////////////////////
void ManyTinyBubbles::generateMoreBubblesFromEmitter()
{
	// TODO: create a standard means to add new bubbles to m_bubbles b/c that operation is done in this method as well as in breakupBubbles()
	// TODO: make sure we aren't recomputing emission positions unnecessarily
	// TODO: consider case where user starts a simulation and then moves the emitter mesh

	//if ( EMITTER_MELTING_RATE != 0 || m_emitter.getSourcePosListSize() == 0 ) {
	//	m_emitter.createEmissionPositionsOnMesh( m_emitter.getLevelSetRes(), m_emitter.getMeltingRate() );
	//}
	m_emitter.createEmissionPositionsOnMesh( m_emitter.getLevelSetRes(), m_emitter.getMeltingRate() );

	// generate bubbles
	std::vector<vec3> new_bubble_positions;
	std::vector<vec3> new_bubble_velocities;
	std::vector<unsigned int> new_bubble_radius_group;
	m_emitter.generateBubbles( m_bubbles.getRadiiList(),
							   new_bubble_positions,
							   new_bubble_velocities,
							   new_bubble_radius_group );

	// add new bubble positions and velocities to m_bubbles
	m_bubbles.addBubblePosToRadiusGroupAtIndex( new_bubble_positions, new_bubble_radius_group );
	m_bubbles.addBubbleVelToRadiusGroupAtIndex( new_bubble_velocities, new_bubble_radius_group );
}


////////////////////////////////////////////////////
// updateBubbleVelocities()
////////////////////////////////////////////////////
void ManyTinyBubbles::updateBubbleVelocities()
{
	// update velocity field to match the Maya fluid
	m_fluid_container.updateVelocityField();

	// compute fraction field for fluid based on bubble positions
	// fraction field influences bubble scattering properties which influences their velocities
	computeFractionField();

	// loop through radius groups
	unsigned int num_radius_groups = m_bubbles.getNumRadiusGroups();
	for ( unsigned int radius_group_index = 0; radius_group_index < num_radius_groups; ++radius_group_index ) {

		// loop through bubbles in the current radius group
		unsigned int num_bubbles_in_list = m_bubbles.getNumBubblesInListWithIndex( radius_group_index );
		for ( unsigned int bubble_index = 0; bubble_index < num_bubbles_in_list; ++bubble_index ) {

			vec3 bubble_pos = m_bubbles.getPosAtIndex( radius_group_index, bubble_index );

			// TODO: tune this b/c it's kind of random

			vec3 bubble_vel_at_prev_time_step = m_bubbles.getVelocityAtIndex( radius_group_index, bubble_index );
			vec3 bubble_vel = vec3( 0.0,
									bubble_vel_at_prev_time_step.Length(),
									0.0 );

			double scattering_probability = computeScatteringProbabilityOfBubble( bubble_vel, bubble_pos );

			// bubble will scatter, so alter bubble direction
			if ( scattering_probability > Convenience::generateRandomDoubleBetweenZeroAndOneInclusive() ) {
				double altered_angle = computeAlteredAngle();

				// altered angle is not zero, so update bubble velocity direction
 				if ( altered_angle < -DBL_EPSILON || altered_angle > DBL_EPSILON ) {

					// update bubble velocity due to scattering
					bubble_vel = computeBubbleVelAfterScattering( bubble_vel, altered_angle );
				}
			}

			// TODO: tune this b/c it's kind of random

			vec3 fluid_vel_at_voxel = m_fluid_container.getVelocityOfVoxelAtPos( bubble_pos );
			vec3 new_bubble_vel = bubble_vel + fluid_vel_at_voxel + vec3( 0.0, 1.0, 0.0 );

			// check for "terminal" velocity
			if ( new_bubble_vel.Length() > 5.0 ) {
				new_bubble_vel /= new_bubble_vel.Length() * 5.0;
			}

			m_bubbles.setVelocityAtIndex( new_bubble_vel,
										  radius_group_index,
										  bubble_index );
		}
	}
}


////////////////////////////////////////////////////
// computeFractionField()
////////////////////////////////////////////////////
void ManyTinyBubbles::computeFractionField()
{
	// set fraction field value at each voxel to 1.0f
	m_fluid_container.resetFractionField();

	// get positions for every bubble present in the simulation
	//std::vector<std::vector<vec3>> bubble_pos_list = m_bubbles.getPosList();

	// loop through radius groups
	unsigned int num_radius_groups = m_bubbles.getNumRadiusGroups();
	for ( unsigned int radius_group_index = 0; radius_group_index < num_radius_groups; ++radius_group_index ) {

		// loop through bubbles in the current radius group
		unsigned int num_bubbles_in_list = m_bubbles.getNumBubblesInListWithIndex( radius_group_index );
		for ( unsigned int bubble_index = 0; bubble_index < num_bubbles_in_list; ++bubble_index ) {

			// reduce fraction field of the voxel the current bubble is inside
			vec3 bubble_pos = m_bubbles.getPosAtIndex( radius_group_index, bubble_index );

			m_fluid_container.reduceFractionFieldOfVoxelAtPos( bubble_pos, m_bubbles.getRadiusAtIndex( radius_group_index ) );
		}
	}

	//// iterate through the bubble radius groups
	//unsigned int radius_index = 0;
	//for ( std::vector<std::vector<vec3>>::iterator outer_it = bubble_pos_list.begin() ; outer_it != bubble_pos_list.end(); ++outer_it ) {
	//	std::vector<vec3> bubble_pos_sublist = *outer_it;

	//	// iterate through the list of vec3s in the current radius group
	//	for ( std::vector<vec3>::iterator inner_it = bubble_pos_sublist.begin() ; inner_it != bubble_pos_sublist.end(); ++inner_it ) {
	//		vec3 bubble_pos = *inner_it;

	//		// reduce fraction field of the voxel the current bubble is inside
	//		m_fluid_container.reduceFractionFieldOfVoxelAtPos( bubble_pos, m_bubbles.getRadiusAtIndex( radius_index ) );
	//	}

	//	++radius_index;
	//}
}


////////////////////////////////////////////////////
// compute scattering probability function, s(x) in paper, [0, 1]
////////////////////////////////////////////////////
double ManyTinyBubbles::computeScatteringProbabilityOfBubble( const vec3& bubble_vel, const vec3& bubble_pos ) const
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
vec3 ManyTinyBubbles::computeBubbleVelAfterScattering( const vec3& old_vel, const double& altered_dir ) const
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
	if ( old_vel[VZ] < -DBL_EPSILON || old_vel[VZ] > DBL_EPSILON ) {

		// generate random numbers for x and y components [-1, 1] so that they are not both 0
		do {
			rotation_axis[VX] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
			rotation_axis[VY] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
		} while ( ( rotation_axis[VX] > -DBL_EPSILON && rotation_axis[VX] < DBL_EPSILON ) &&
				  ( rotation_axis[VY] > -DBL_EPSILON && rotation_axis[VY] < DBL_EPSILON ) );

		rotation_axis[VZ] = -1.0 * ( rotation_axis[VX] * old_vel[VX] + rotation_axis[VY] * old_vel[VY] ) / old_vel[VZ];
	}
	else if ( old_vel[VY] < -DBL_EPSILON || old_vel[VY] > DBL_EPSILON ) {

		// generate random numbers for x and z components [-1, 1] so that they are not both 0
		do {
			rotation_axis[VX] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
			rotation_axis[VZ] = Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive();
		} while ( ( rotation_axis[VX] > -DBL_EPSILON && rotation_axis[VX] < DBL_EPSILON ) &&
				  ( rotation_axis[VZ] > -DBL_EPSILON && rotation_axis[VZ] < DBL_EPSILON ) );

		rotation_axis[VZ] = -1.0 * ( rotation_axis[VX] * old_vel[VX] + rotation_axis[VZ] * old_vel[VZ] ) / old_vel[VY];
	}
	else if ( old_vel[VX] < -DBL_EPSILON || old_vel[VX] > DBL_EPSILON ) {

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
// breakupBubbles()
////////////////////////////////////////////////////
void ManyTinyBubbles::breakupBubbles()
{
	// TODO: create a standard means to add new bubbles to m_bubbles b/c that operation is done in this method as well as in generateMoreBubblesFromEmitter()

	std::vector<BubbleLocator> bubbles_to_remove;
	std::vector<BubbleStruct> bubbles_to_add;

	// loop through radius groups
	// starting index is 1 b/c bubble cannot split if it belongs to the smallest radius group
	unsigned int num_radius_groups = m_bubbles.getNumRadiusGroups();
	for ( unsigned int radius_group_index = 1; radius_group_index < num_radius_groups; ++radius_group_index ) {

		// loop through bubbles in the current radius group
		unsigned int num_bubbles_in_list = m_bubbles.getNumBubblesInListWithIndex( radius_group_index );
		for ( unsigned int bubble_index = 0; bubble_index < num_bubbles_in_list; ++bubble_index ) {

			// bubble wants to split
			if ( m_bubbles.getBreakupFrequency() > Convenience::generateRandomDoubleBetweenZeroAndOneInclusive() ) {

				// get bubble's current position and velocity
				vec3 bubble_pos = m_bubbles.getPosAtIndex( radius_group_index, bubble_index );
				vec3 bubble_vel = m_bubbles.getVelAtIndex( radius_group_index, bubble_index );

				// compute new positions for two new bubbles created from split bubble
				vec3 new_pos_1, new_pos_2;
				splitBubble( bubble_pos,
							 m_bubbles.getRadiusAtIndex( radius_group_index ),
							 new_pos_1,
							 new_pos_2 );

				bubbles_to_add.push_back( BubbleStruct( new_pos_1,
														bubble_vel,
														radius_group_index - 1 ) );
				bubbles_to_add.push_back( BubbleStruct( new_pos_2,
														bubble_vel,
														radius_group_index - 1 ) );

				// mark bubble that just split to be removed from bubble list
				bubbles_to_remove.push_back( BubbleLocator( radius_group_index, bubble_index ) );
			}
		}
	}

	// add new bubbles m_bubbles data structure
	for ( std::vector<BubbleStruct>::iterator it = bubbles_to_add.begin(); it != bubbles_to_add.end(); ++it ) {
		m_bubbles.addBubble( it->radius_group,
							 it->pos,
							 it->vel );
	}

	deleteBubblesInList( bubbles_to_remove );
}


////////////////////////////////////////////////////
// compute new positions for two new bubbles created from split bubble
////////////////////////////////////////////////////
void ManyTinyBubbles::splitBubble( const vec3&		current_pos,
								   const double&	current_radius,
								   vec3&			new_pos_1,
								   vec3&			new_pos_2 ) const
{
	// TODO: remove "polySurface1" hardcoded value
	// TODO: move hardcoded "1000" value to top of this file

	// random angle for use when generating initial positions of new bubbles formed from split
	double phi = Convenience::generateRandomDoubleInclusive( 0.0, 2.0 ) * M_PI;

	// generate new positions located at the two end points of horizontal diameter across parent bubble
	new_pos_1 = current_pos + vec3( current_radius * cos(phi),
									0.0,
									current_radius * sin( phi ) );
	new_pos_2 = current_pos - vec3( current_radius * cos( phi ),
									0.0,
									current_radius * sin( phi ) );

	// check if newly generated positions are outside fluid container
	// if they are, simply assign them the position of the bubble before it split

	//if ( m_fluid_container.posIsOutsideFluidContainer( new_pos_1 ) ) {
	//	new_pos_1 = current_pos;
	//}
	//if ( m_fluid_container.posIsOutsideFluidContainer( new_pos_2 ) ) {
	//	new_pos_2 = current_pos;
	//}

	int exists;
	MGlobal::executeCommand( "objExists polySurface1", exists );

	// if polySurface1 exists, then fluid has been converted into a polygon mesh so we use the level set method to determine fluid surface boundary
	if ( exists == 1 ) {

		// if too many bubbles, then filter bubbles so we only concern ourselves with bubbles near the fluid's surface
		if ( m_bubbles.getTotalBubbleNumber() > 1000 ) {

			////////////////////////////////////////////////////
			// compute for pos_1
			////////////////////////////////////////////////////

			// bubbles are outside fluid surface boundary
			if ( m_fluid_container.posIsUnderBoundary( new_pos_1 ) == 1 ) {
				new_pos_1 = current_pos;
			}
			// bubbles are located near the fluid surface boundary
			else if ( m_fluid_container.posIsUnderBoundary( new_pos_1 ) == 0 ) {

				// use finer level set method to determine if bubbles should be removed
				if ( m_fluid_container.isPosOutsideFluidViaLevelSet( new_pos_1 ) ) {
					new_pos_1 = current_pos;
				}
			}


			////////////////////////////////////////////////////
			// compute for pos_2
			////////////////////////////////////////////////////

			if ( m_fluid_container.posIsUnderBoundary( new_pos_2 ) == 1 ) {
				new_pos_2 = current_pos;
			}
			else if ( m_fluid_container.posIsUnderBoundary( new_pos_2 ) == 0 ) {
				if ( m_fluid_container.isPosOutsideFluidViaLevelSet( new_pos_2 ) ) {
					new_pos_2 = current_pos;
				}
			}
		}
		// if there aren't too many bubbles, directly use the finer level set method to filter bubbles
		else {

			// compute for pos_1
			if ( m_fluid_container.isPosOutsideFluidViaLevelSet( new_pos_1 ) ) {
				new_pos_1 = current_pos;
			}

			// compute for pos_2
			if ( m_fluid_container.isPosOutsideFluidViaLevelSet( new_pos_2 ) ) {
				new_pos_2 = current_pos;
			}
		}
	}
	else {

		// compute for pos_1
		if ( m_fluid_container.isPosOutsideFluidViaDensity( new_pos_1 ) ) {
			new_pos_1 = current_pos;
		}

		// compute for pos_2
		if ( m_fluid_container.isPosOutsideFluidViaDensity( new_pos_2 ) ) {
			new_pos_2 = current_pos;
		}
	}
}


////////////////////////////////////////////////////
// store the first emitter mesh information (vertext list, face list) in GlobalState
////////////////////////////////////////////////////
void ManyTinyBubbles::storeMeshInfoByName( MString emitter_mesh_name ) const
{
	// select Maya object by name
	MSelectionList maya_sel_list;
	MGlobal::getSelectionListByName( emitter_mesh_name, maya_sel_list );

	// get path to a mesh DAG node
	MDagPath dag_node_path;
	maya_sel_list.getDagPath( 0, dag_node_path );

	// get MFnMesh from MDagPath
	MFnMesh mesh_surface( dag_node_path );

	// get number of triangles in mesh and the vertex indices for each triangle
	// triangle_vertex_indices will be three times larger than triangle_num
	MIntArray minta_tri_num;
	MIntArray minta_tri_vertex_indices;
	mesh_surface.getTriangles( minta_tri_num,
							   minta_tri_vertex_indices );
	unsigned int triangle_vertex_num = minta_tri_vertex_indices.length();

	// copy triangle_vertex_indices into a C++ vector
	std::vector<int> triangle_vertex_indices;
	triangle_vertex_indices.resize( triangle_vertex_num );
	minta_tri_vertex_indices.get( &triangle_vertex_indices[0] );


	////////////////////////////////////////////////////
	// create face list
	////////////////////////////////////////////////////

	std::vector<int> face_list;
	for ( unsigned int i = 0; i < triangle_vertex_num; i += 3 ) {
		face_list.push_back( triangle_vertex_indices[i] + 1 );
		face_list.push_back( triangle_vertex_indices[i + 1] + 1 );
		face_list.push_back( triangle_vertex_indices[i + 2] + 1 );
	}

	// store face list data in GlobalState
	GlobalState::storeFaceList( face_list );
	

	////////////////////////////////////////////////////
	// create vertex list
	////////////////////////////////////////////////////

	std::vector<vec3> vert_list;

	MFloatPointArray mesh_vertices;
	mesh_surface.getPoints(mesh_vertices, MSpace::kWorld);

	for ( unsigned int i = 0; i < mesh_vertices.length(); ++i ) {
		MFloatPoint vertex = mesh_vertices[i];
		vec3 new_vert( vertex[VX],
						vertex[VY],
						vertex[VZ] );
		vert_list.push_back( new_vert );
	}

	// store vertex list data in GlobalState
	GlobalState::storePointList( vert_list );
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

	//ManyTinyBubbles::m_output = nAttr.create( "output", "out", MFnNumericData::kFloat, 0.0 );
	//nAttr.setWritable( true );
	//nAttr.setStorable( false );
	//nAttr.setKeyable( false );

	ManyTinyBubbles::m_output = tAttr.create(  "output_mesh", "out", MFnData::kMesh, &stat );
	nAttr.setWritable( true );
	nAttr.setStorable( false );
	nAttr.setKeyable( false );

	ManyTinyBubbles::m_time = uAttr.create( "time", "tm", MFnUnitAttribute::kTime, 0.0 );
	uAttr.setWritable( true );
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

	ManyTinyBubbles::m_fluid_transform_name = tAttr.create( "fluid_transform_name", "ftn", MFnData::kString );
 	tAttr.setStorable( true );
    tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::m_node_name = tAttr.create( "node_name", "nn", MFnData::kString );
 	tAttr.setStorable( true );
    tAttr.setWritable( false );
	tAttr.setStorable( false );

	ManyTinyBubbles::m_melting_rate = nAttr.create( "melting_rate", "mr", MFnNumericData::kInt, 0 );
	nAttr.setMin( 0 );
	nAttr.setMax( 100 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );

	ManyTinyBubbles::m_level_set_resolution = nAttr.create( "level_set_resolution", "lsr", MFnNumericData::kInt, 1000 );
	nAttr.setMin( 100 );
	nAttr.setMax( 10000 );
 	nAttr.setStorable( true );
	nAttr.setWritable( true );
 	nAttr.setKeyable( true );


	////////////////////////////////////////////////////
	// add attributes to node
	////////////////////////////////////////////////////

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

	stat = addAttribute( ManyTinyBubbles::m_fluid_transform_name );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_node_name );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_melting_rate );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }
	stat = addAttribute( ManyTinyBubbles::m_level_set_resolution );
		if ( !stat ) { stat.perror( "addAttribute" ); return stat; }


	////////////////////////////////////////////////////
	// setup dependencies between attributes
	////////////////////////////////////////////////////

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

	stat = attributeAffects( ManyTinyBubbles::m_fluid_transform_name, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_node_name, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_melting_rate, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }
	stat = attributeAffects( ManyTinyBubbles::m_level_set_resolution, ManyTinyBubbles::m_output );
		if ( !stat ) { stat.perror( "attributeAffects" ); return stat; }

	return MS::kSuccess;
}