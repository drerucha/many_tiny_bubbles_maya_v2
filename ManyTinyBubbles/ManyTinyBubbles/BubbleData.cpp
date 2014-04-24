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
}


////////////////////////////////////////////////////
// setRadii()
////////////////////////////////////////////////////
void BubbleData::setRadii( const double& radius_min,
						   const double& radius_max )
{
	// TODO: add back logic that ensures bubble radii are generated so each bubble size it twice the volume of the size immediately smaller than it

	// clear m_radii_list
	if ( m_radii_list.size() > 0 ) {
		m_radii_list.clear();
	}

	double diff = abs( radius_max - radius_min );
	double step = diff / ( NUM_RADII - 1 );

	// fill m_radii_list
	for ( unsigned int i = 0; i < NUM_RADII; ++i ) {
		m_radii_list.push_back( radius_min + step * i );
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
	m_radii_list.clear();
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
	std::vector<vec3> pos_list = m_pos_list[i];
	m_pos_list.at( i ).erase( pos_list.begin() + j );
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

//std::vector<vec3> BubbleData::getPosListForRadiusGroupAtIndex( const unsigned int& i ) const
//{
//	return m_pos_list.at( i );
//}