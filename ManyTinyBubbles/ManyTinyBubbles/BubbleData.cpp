#include "BubbleData.h"

#include <maya/MGlobal.h>

#include "Convenience.h"



////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

BubbleData::BubbleData()
{
}

BubbleData::~BubbleData()
{
}


////////////////////////////////////////////////////
// initialize member variables
////////////////////////////////////////////////////
void BubbleData::init( double scattering_frequency,
					   double scattering_coefficient,
					   double breakup_frequency,
					   double size_min,
					   double size_max )
{
	m_scattering_frequency = scattering_frequency;
	m_scattering_coefficient = scattering_coefficient;
	m_breakup_frequency = breakup_frequency;
	m_size_min = size_min;
	m_size_max = size_max;

	// populate m_radii_list
	setRadii( size_min, size_max );


	// TODO: don't repeat this block of code here and in reset()

	// TODO: this might be a dumb thing to do
	// init m_pos_list and m_vel_list with empty std::vector<vec3>
	// only perform this action once
	if ( m_pos_list.size() == 0 ) {
		for ( unsigned int i = 0; i < m_radii_list.size(); ++i ) {
			std::vector<vec3> empty_pos, empty_vel;
			m_pos_list.push_back( empty_pos );
			m_vel_list.push_back( empty_vel );
		}
	}
}


////////////////////////////////////////////////////
// setRadii()
////////////////////////////////////////////////////
void BubbleData::setRadii( const double& radius_min,
						   const double& radius_max )
{
	// clear m_radii_list
	if ( m_radii_list.size() > 0 ) {
		m_radii_list.clear();
	}

	//double diff = abs( radius_max - radius_min );
	//double step = diff / ( NUM_RADII - 1 );

	//// fill m_radii_list
	//for ( unsigned int i = 0; i < NUM_RADII; ++i ) {
	//	m_radii_list.push_back( radius_min + step * i );
	//}


	// generate possible bubble radii in such a way that each size up is twice the volume of the previous size

	double l_radius_max = radius_max;
	double l_radius_min = radius_min;

	// ensure radius_max > radius_min
	if ( l_radius_max < l_radius_min ) {
		double tmp = l_radius_min;
		l_radius_min = l_radius_max;
		l_radius_max = tmp;
	}
	
	unsigned int segment = 1;
	while ( l_radius_max / l_radius_min >= 2.0 ) {
		l_radius_max = l_radius_max / 2.0;
		++segment;
	}
	
	const double BASE = pow( 2.0, 1.0 / 3.0 );
	for ( unsigned int i = 0; i < segment; ++i ) {
		double radius = l_radius_min * pow( BASE, ( double )i );
		m_radii_list.push_back( radius );
	}
}


////////////////////////////////////////////////////
// deleteAllParticlesInMaya()
////////////////////////////////////////////////////
void BubbleData::deleteAllParticlesInMaya()
{
	// particles are separated into groups based on their radii
	for ( unsigned int i = 1; i <= getNumRadii(); ++i ) {
		int particle_exists = checkIfParticleExists( i );
		if ( particle_exists ) {
			deleteParticle( i );
		}
	}
}


////////////////////////////////////////////////////
// check if particle exists in Maya
////////////////////////////////////////////////////
int BubbleData::checkIfParticleExists( const unsigned int& num ) const
{
	// TODO: change "bubbleParticle" to be the constant string set in ManyTinyBubblesNode
	// maybe make the name a member variable of BubbleData

	std::string cmd = "particleExists bubbleParticle";
	Convenience::appendNumToStdString( cmd, num );
	cmd += ";";

	int particle_exists;
	MGlobal::executeCommand( Convenience::convertStdStringToMString( cmd ), particle_exists );

	return particle_exists;
}


////////////////////////////////////////////////////
// select particle then delete it in Maya
////////////////////////////////////////////////////
void BubbleData::deleteParticle( const unsigned int& num ) const
{
	// TODO: change "bubbleParticle" to be the constant string set in ManyTinyBubblesNode
	// maybe make the name a member variable of BubbleData

	std::string cmd = "select -replace bubbleParticle";
	Convenience::appendNumToStdString( cmd, num );
	cmd += ";";

	MGlobal::executeCommand( Convenience::convertStdStringToMString( cmd ) );
	MGlobal::executeCommand( "doDelete" );
}


////////////////////////////////////////////////////
// reset
////////////////////////////////////////////////////
void BubbleData::reset()
{
	m_pos_list.clear();
	m_vel_list.clear();
		 
	for ( unsigned int i = 0; i < m_radii_list.size(); ++i ) {
		std::vector<vec3> empty_pos, empty_vel;
		m_pos_list.push_back( empty_pos );
		m_vel_list.push_back( empty_vel );
	}
}


////////////////////////////////////////////////////
// i indicates radius group
// j indicates position within a specific radius group
////////////////////////////////////////////////////
void BubbleData::removeBubbleAtIndex( const unsigned int& i,
									  const unsigned int& j )
{
	// TODO: test this

	// m_pos_list is std::vector<std::vector<vec3>>
	//std::vector<vec3> pos_list = m_pos_list[i];
	//m_pos_list.at( i ).erase( pos_list.begin() + j );

	// m_vel_list is std::vector<std::vector<vec3>>
	//std::vector<vec3> vel_list = m_vel_list[i];
	//m_vel_list.at( i ).erase( vel_list.begin() + j );

	// m_pos_list is std::vector<std::vector<vec3>>
	// m_vel_list is std::vector<std::vector<vec3>>
	m_pos_list[i].erase( m_pos_list[i].begin() + j );
	m_vel_list[i].erase( m_vel_list[i].begin() + j );
}


////////////////////////////////////////////////////
// i indicates radius group
////////////////////////////////////////////////////
void BubbleData::addBubblePosToRadiusGroupAtIndex( const vec3& pos,
												   const unsigned int& i )
{
	// TODO: test this

	// m_pos_list is std::vector<std::vector<vec3>>
	m_pos_list.at( i ).push_back( pos );
}


////////////////////////////////////////////////////
// i indicates radius group
////////////////////////////////////////////////////
void BubbleData::addBubblePosToRadiusGroupAtIndex( std::vector<vec3>			pos_list,
												   std::vector<unsigned int>	radius_group_index_list )
{
	// TODO: test this

	if ( pos_list.size() != radius_group_index_list.size() ) {
		Convenience::printInScriptEditor( MString( "ERROR: vectors are not same size in BubbleData::addBubblePosToRadiusGroupAtIndex" ) );
	}
	else {
		for ( unsigned int i = 0; i < pos_list.size(); ++i ) {
			unsigned int radius_group_index = radius_group_index_list[i];
			m_pos_list[radius_group_index].push_back( pos_list[i] );
		}
	}
}


////////////////////////////////////////////////////
// i indicates radius group
////////////////////////////////////////////////////
void BubbleData::addBubbleVelToRadiusGroupAtIndex( std::vector<vec3>			vel_list,
												   std::vector<unsigned int>	radius_group_index_list )
{
	// TODO: test this

	if ( vel_list.size() != radius_group_index_list.size() ) {
		Convenience::printInScriptEditor( MString( "ERROR: vectors are not same size in BubbleData::addBubbleVelToRadiusGroupAtIndex" ) );
	}
	else {
		for ( unsigned int i = 0; i < vel_list.size(); ++i ) {
			unsigned int radius_group_index = radius_group_index_list[i];
			m_vel_list[radius_group_index].push_back( vel_list[i] );
		}
	}
}


////////////////////////////////////////////////////
// create Maya particle groups from m_pos_list
////////////////////////////////////////////////////
void BubbleData::createMayaParticlesWithName( const std::string& particle_name ) const
{
	for ( unsigned int radius_group_index = 0; radius_group_index < m_pos_list.size(); ++radius_group_index ) {
		std::vector<vec3> bubble_pos_list = m_pos_list.at( radius_group_index );

		if ( bubble_pos_list.size() > 0 ) {
			std::string create_cmd = "particle ";

			for ( std::vector<vec3>::iterator it = bubble_pos_list.begin(); it != bubble_pos_list.end(); ++it ) {
				create_cmd += "-position ";
				vec3 pos = *it;
				Convenience::appendVec3ToStdString( create_cmd, pos );
				create_cmd += " ";
			}

			create_cmd += "-conserve 1 -name ";
			create_cmd += particle_name;
			Convenience::appendNumToStdString( create_cmd, radius_group_index + 1 );

			// create particle
			// "particle [-position x y z]+ -conserve 1 -name bubbleParticle1"
			MGlobal::executeCommand( Convenience::convertStdStringToMString( create_cmd ) );

			// set particle to render as sphere
			std::string object_name = particle_name;
			Convenience::appendNumToStdString( object_name, radius_group_index + 1 );
			Convenience::setParticleRenderTypeToSphere( Convenience::convertStdStringToMString( object_name ) );
		}
	}
}


////////////////////////////////////////////////////
// set radii of Maya particles
////////////////////////////////////////////////////
void BubbleData::setRadiiForMayaParticlesWithName( const std::string& particle_name ) const
{
	for ( unsigned int i = 0; i < m_radii_list.size(); ++i ) {
		int particle_exists = checkIfParticleExists( i + 1 );
		if ( particle_exists ) {
			std::string object_name = particle_name;
			Convenience::appendNumToStdString( object_name, i + 1 );

			// TODO: ask about this long MEL command
			// addAttr -internalSet true -longName radius -attributeType "float" -minValue 0 -maxValue 10 -defaultValue 0.5 bubbleParticle1Shape
			std::string cmd = "addAttr -internalSet true -longName radius -attributeType \"float\" -minValue 0 -maxValue 10 -defaultValue 0.5 ";
			cmd += object_name;
			cmd += "Shape";
			MGlobal::executeCommand( Convenience::convertStdStringToMString( cmd ) );

			MGlobal::executeCommand( "setAttr " + Convenience::convertStdStringToMString( object_name ) + "Shape.radius " + m_radii_list[i] );
		}
	}
}


////////////////////////////////////////////////////
// set bubble velocity at specific index
// i indicates radius group
// j indicates position within a specific radius group
////////////////////////////////////////////////////
void BubbleData::setVelocityAtIndex( const vec3& new_vel,
									 const unsigned int& i,
									 const unsigned int& j )
{
	// m_vel_list is std::vector<std::vector<vec3>>
	m_vel_list[i][j] = new_vel;
}


////////////////////////////////////////////////////
// update bubble position at index using corresponding velocity
////////////////////////////////////////////////////
//void BubbleData::updatePosAtIndex( const unsigned int& i, const unsigned int& j )
//{
//	m_pos_list[i][j] += 
//}


////////////////////////////////////////////////////
// performs explicit Euler integration to update all bubble positions using their velocities and the time step
////////////////////////////////////////////////////
void BubbleData::updateBubblePositions( const float& time_step )
{
	for ( unsigned int i = 0; i < m_pos_list.size(); ++i ) {
		std::vector<vec3> sublist = m_pos_list[i];

		for ( unsigned int j = 0; j < sublist.size(); ++j ) {
			m_pos_list[i][j] += time_step * m_vel_list[i][j];
		}
	}
}


////////////////////////////////////////////////////
// performs explicit Euler integration to update a single bubble's position using its velocity and the time step
////////////////////////////////////////////////////
void BubbleData::updateBubblePositionsAtIndex( const unsigned int&	i,
											   const unsigned int&	j,
											   const float&			time_step )
{
	m_pos_list[i][j] += time_step * m_vel_list[i][j];
}


////////////////////////////////////////////////////
// add new bubble
// add new vec3 position to m_pos_list and vec3 velocity to m_vel_list
////////////////////////////////////////////////////
void BubbleData::addBubble( const unsigned int&	radius_group_index,
							vec3				pos,
							vec3				vel )
{
	if ( radius_group_index > m_radii_list.size() ) {
		Convenience::printInScriptEditor( MString( "ERROR: radius_group_index greater than number of radius groups in BubbleData::addBubbleToRadiusGroup" ) );
	}
	else {
		m_pos_list[radius_group_index].push_back( pos );
		m_vel_list[radius_group_index].push_back( vel );
	}
}


////////////////////////////////////////////////////
// getters
////////////////////////////////////////////////////

unsigned int BubbleData::getNumRadii() const
{
	return ( unsigned int )m_radii_list.size();
}

std::vector<std::vector<vec3>> BubbleData::getPosList() const
{
	return m_pos_list;
}

double BubbleData::getRadiusAtIndex( const unsigned int& index ) const
{
	return m_radii_list[index];
}

double BubbleData::getScatteringFrequency() const
{
	return m_scattering_frequency;
}

double BubbleData::getScatteringCoefficient() const
{
	return m_scattering_coefficient;
}

double BubbleData::getBreakupFrequency() const
{
	return m_breakup_frequency;
}

std::vector<double> BubbleData::getRadiiList() const
{
	return m_radii_list;
}

vec3 BubbleData::getVelocityAtIndex( const unsigned int& i, const unsigned int& j ) const
{
	return m_vel_list[i][j];
}

//std::vector<vec3> BubbleData::getPosListForRadiusGroupAtIndex( const unsigned int& i ) const
//{
//	return m_pos_list.at( i );
//}

unsigned int BubbleData::getNumRadiusGroups() const
{
	return ( unsigned int )m_radii_list.size();
}

unsigned int BubbleData::getNumBubblesInListWithIndex( const unsigned int& i ) const
{
	return ( unsigned int ) m_pos_list[i].size();
}

vec3 BubbleData::getPosAtIndex( const unsigned int& i, const unsigned int& j ) const
{
	return m_pos_list[i][j];
}

vec3 BubbleData::getVelAtIndex( const unsigned int& i, const unsigned int& j ) const
{
	return m_vel_list[i][j];
}

unsigned int BubbleData::getTotalBubbleNumber() const
{
	unsigned int total_bubble_num = 0;
    for( unsigned int i = 0; i < m_radii_list.size(); ++i ) {
		total_bubble_num +=  ( unsigned int )m_pos_list[i].size();
	}
	return total_bubble_num;
}