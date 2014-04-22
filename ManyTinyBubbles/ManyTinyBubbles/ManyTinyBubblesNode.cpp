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

		// TODO: set fraction field to some default value, 1.0f
		// TODO: allow fraction field to handle voxels that aren't squares
		//mFractionField.initialize( 1.0f, m_fluid_container_cell_size_x, m_fluid_container_dim_x, m_fluid_container_dim_y, m_fluid_container_dim_z );

		// TODO: compute fraction field
		//computeFractionField();

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
// simulate bubbles
////////////////////////////////////////////////////
void ManyTinyBubbles::advectParticles( const float& dt )
{
	// update velocity field to match the Maya fluid
	m_fluid_container.updateVelocityField();

	// get positions for every bubble present in the simulation
	std::vector<std::vector<vec3>> bubble_pos_list = m_bubbles.getPosList();

	// iterate through the list of bubble position lists (each bubble radius has a unique list of positions)
	for ( std::vector<std::vector<vec3>>::iterator outer_it = bubble_pos_list.begin() ; outer_it != bubble_pos_list.end(); ++outer_it ) {
		std::vector<vec3> bubble_pos_sublist = *outer_it;

		// iterate through the list of vec3s in one of the radius groups (one of the sublists of bubble_pos_list)
		for ( std::vector<vec3>::iterator inner_it = bubble_pos_sublist.begin() ; inner_it != bubble_pos_sublist.end(); ++inner_it ) {
			vec3 bubble_pos = *inner_it;

			// get velocity of cell in fluid container
			vec3 cell_vel = m_fluid_container.getVelocityAtPos( bubble_pos );


		}
	}





	//int j = 0;
	//for ( std::vector<std::vector<vec3>>::iterator iterRadius = bubblePosList.begin(); iterRadius != bubblePosList.end(); ++iterRadius, ++j ) {
	//	int i = 0;
	//	for ( std::vector<vec3>::iterator iterPos = bubblePosList[j].begin(), iterVel = bubbleVelList[j].begin(); iterPos != bubblePosList[j].end(); ++i ) {
	//		vec3 position = bubblePosList[j][i];

	//		int position_grid_X = (int)((position[0] + CONTAINER_DIM_X * CELL_SIZE / 2.0f - CONTAINER_TRANS_X) / CELL_SIZE);
	//		int position_grid_Y = (int)((position[1] + CONTAINER_DIM_Y * CELL_SIZE / 2.0f - CONTAINER_TRANS_Y) / CELL_SIZE);
	//		int position_grid_Z = (int)((position[2] + CONTAINER_DIM_Z * CELL_SIZE / 2.0f - CONTAINER_TRANS_Z) / CELL_SIZE);



	//		vec3 velocity(velocityArray[position_grid_X + position_grid_Y * CONTAINER_DIM_X + position_grid_Z * CONTAINER_DIM_X * CONTAINER_DIM_Y + 0],
	//			          velocityArray[position_grid_X + position_grid_Y * CONTAINER_DIM_X + position_grid_Z * CONTAINER_DIM_X * CONTAINER_DIM_Y + 1],
	//					  velocityArray[position_grid_X + position_grid_Y * CONTAINER_DIM_X + position_grid_Z * CONTAINER_DIM_X * CONTAINER_DIM_Y + 2]);

	//		double ran_num = ( rand() % 100 ) / 100.0f;
	//		double velocityMag = velocity.Length(); 
	//		double fractionField = mFractionField( position_grid_X, position_grid_Y, position_grid_Z );
	//		double scatterOdd = scatterFreq * 100* ( 1 - fractionField ) * velocityMag * velocityMag;	//between 0~1 

	//		// alter the direction
	//		if ( scatterOdd > ran_num ) {
	//			double x = 2 * ran_num * scatterCoef - scatterCoef + 1;
	//			double y = 2 * ran_num + scatterCoef - 1;
	//			double cosTheta;
	//			double theta;

	//			if ( x == 0 ) {
	//				theta = PI / 2.0f;
	//			}
	//			else {
	//				cosTheta = y / x;
	//				theta = acos( cosTheta );
	//			}

 //				if ( theta != 0 ) {
	//				double rotateAxisX = 1;
	//				double rotateAxisY = 0;
	//				double rotateAxisZ = 0;

	//				if ( velocity[2] != 0 ) {
	//					rotateAxisX = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//					rotateAxisY = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;

	//					while ( rotateAxisX == 0 && rotateAxisY == 0 ) {
	//						rotateAxisX = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//						rotateAxisY = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//					}

	//					rotateAxisZ = -( rotateAxisX * velocity[0] + rotateAxisY * velocity[1] ) / velocity[2];
	//				}
	//				else if ( velocity[1] != 0 ) {
	//					rotateAxisX = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//					rotateAxisZ = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;

	//					while ( rotateAxisX == 0 && rotateAxisZ == 0 ) {
	//						rotateAxisX = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//						rotateAxisZ = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//					}

	//					rotateAxisY = -( rotateAxisX * velocity[0] + rotateAxisZ * velocity[2] ) / velocity[1];
	//				}
	//				else if ( velocity[0] != 0 ) {
	//					rotateAxisY = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//					rotateAxisZ = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;

	//					while ( rotateAxisY == 0 && rotateAxisZ == 0 ) {
	//						rotateAxisY = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//						rotateAxisZ = ( ( ( rand() % 200 ) + 1 ) - 100 ) / 100.0f;
	//					}

	//					rotateAxisX = -( rotateAxisY * velocity[1] + rotateAxisZ * velocity[2] ) / velocity[0];
	//				}

	//				double length = sqrt( rotateAxisX * rotateAxisX + rotateAxisY * rotateAxisY + rotateAxisZ * rotateAxisZ );
	//				rotateAxisX = rotateAxisX / length;
	//				rotateAxisY = rotateAxisY / length;
	//				rotateAxisZ = rotateAxisZ / length;

	//				double ss = cos( theta / 2.0f );
	//				double xx = sin( theta / 2.0f ) * rotateAxisX;
	//				double yy = sin( theta / 2.0f ) * rotateAxisY;
	//				double zz = sin( theta / 2.0f ) * rotateAxisZ;

	//				double newX = ( 1 - 2 * yy * yy - 2 * zz * zz ) * velocity[0] + ( 2 * xx * yy - 2 * ss * zz ) * velocity[1] + ( 2 * xx * zz + 2 * ss * yy ) * velocity[2];
	//				double newY = ( 2 * xx * yy + 2 * ss * zz ) * velocity[0] + ( 1 - 2 * xx * xx - 2 * zz * zz ) * velocity[1] + ( 2 * yy * zz - 2 * ss * xx ) * velocity[2];
	//				double newZ = ( 2 * xx * zz - 2 * ss * yy ) * velocity[0] + ( 2 * yy * zz + 2 * ss * xx ) * velocity[1] + ( 1 - 2 * xx * xx - 2 * yy * yy ) * velocity[2];
	//				velocity[0] = newX;
	//				velocity[1] = newY;
	//				velocity[2] = newZ;
	//			}

	//		}

	//		position += velocity;

	//		// particle goes outside the container
	//		if ( position[0] - CONTAINER_TRANS_X < -CONTAINER_DIM_X * CELL_SIZE / 2.0f ||
	//			 position[0] - CONTAINER_TRANS_X >  CONTAINER_DIM_X * CELL_SIZE / 2.0f ||
	//			 position[1] - CONTAINER_TRANS_Y < -CONTAINER_DIM_Y * CELL_SIZE / 2.0f ||
	//			 position[1] - CONTAINER_TRANS_Y >= (CONTAINER_DIM_Y * CELL_SIZE) / 2.0f - 0.001 ||
	//			 position[2] - CONTAINER_TRANS_Z < -CONTAINER_DIM_Z * CELL_SIZE / 2.0f ||
	//			 position[2] - CONTAINER_TRANS_Z >  CONTAINER_DIM_Z * CELL_SIZE / 2.0f )
	//		{
	//			if ( iterPos == bubblePosList[j].begin() ) {
	//				bubblePosList[j].erase( iterPos );
	//				iterPos = bubblePosList[j].begin();
	//				bubbleVelList[j].erase( iterVel );
	//				iterVel = bubbleVelList[j].begin();
	//				i--;
	//			}
	//			else {
	//				--( iterPos = bubblePosList[j].erase( iterPos ) );
	//				--( iterVel = bubbleVelList[j].erase( iterVel ) );
	//				i--;
	//				++iterPos; 
	//				++iterVel;
	//			}
	//		}
	//		else {
	//			bubblePosList[j][i] = position;

	//			// bubble break
	//			if ( breakFreq > ( rand() % 100 ) / 100.0f && j > 0 ) {
	//				vec3 newPosition1 = ( position + vec3( bubbleRadiusList[j], 0, 0 ) );
	//				vec3 newPosition2 = ( position - vec3( bubbleRadiusList[j], 0, 0 ) );

	//				if ( newPosition1[0] - CONTAINER_TRANS_X  < -CONTAINER_DIM_X * CELL_SIZE / 2.0f ||
	//					 newPosition1[0] - CONTAINER_TRANS_X  > CONTAINER_DIM_X * CELL_SIZE / 2.0f )
	//				{
	//				    newPosition1 = position;
	//				}

	//				if ( newPosition2[0] - CONTAINER_TRANS_X  < -CONTAINER_DIM_X * CELL_SIZE / 2.0f ||
	//					 newPosition2[0] - CONTAINER_TRANS_X  > CONTAINER_DIM_X * CELL_SIZE / 2.0f )
	//				{	 
	//					 newPosition2 = position;
	//				}

	//				bubblePosList[j-1].push_back( newPosition1 );
	//				bubbleVelList[j-1].push_back( bubbleVelList[j][i] );
	//				bubblePosList[j-1].push_back( newPosition2 );
	//				bubbleVelList[j-1].push_back( bubbleVelList[j][i] );

	//				if ( iterPos == bubblePosList[j].begin() ) {
	//					bubblePosList[j].erase( iterPos );
	//					iterPos = bubblePosList[j].begin();
	//					bubbleVelList[j].erase( iterVel );
	//					iterVel = bubbleVelList[j].begin();
	//					i--;
	//				}
	//				else {
	//					--( iterPos = bubblePosList[j].erase( iterPos ) );
	//					--( iterVel = bubbleVelList[j].erase( iterVel ) );
	//					i--;
	//					++iterPos; 
	//					++iterVel;
	//				}
	//			}
	//			else {
	//				++iterPos; 
	//				++iterVel;
	//			}
	//		}

	//	}
	//}


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