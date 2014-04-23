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