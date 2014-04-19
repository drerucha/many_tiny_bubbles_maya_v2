#include "BubbleData.h"

#include "Convenience.h"
#include <maya/MGlobal.h>


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
void BubbleData::init( float scattering_frequency,
					   float scattering_coefficient,
					   float breakup_frequency,
					   float size_min,
					   float size_max )
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
void BubbleData::setRadii( const float& radius_min,
						   const float& radius_max )
{
	// clear m_radii_list
	if ( m_radii_list.size() > 0 ) {
		m_radii_list.clear();
	}

	float diff = abs( radius_max - radius_min );
	float step = diff / ( NUM_RADII - 1 );

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