#include "BubbleData.h"

#include "Convenience.h"
#include <maya/MGlobal.h>


////////////////////////////////////////////////////
// constructor
////////////////////////////////////////////////////
BubbleData::BubbleData(void)
{
}


////////////////////////////////////////////////////
// destructor
////////////////////////////////////////////////////
BubbleData::~BubbleData(void)
{
}


////////////////////////////////////////////////////
// setRadii()
////////////////////////////////////////////////////
void BubbleData::setRadii( const float& radius_min, const float& radius_max )
{
	m_radii_list.clear();

	float diff = abs( radius_max - radius_min );
	float step = diff / ( NUM_RADII - 1 );

	for ( unsigned int i = 0; i < NUM_RADII; ++i ) {
		m_radii_list.push_back( radius_min + step * i );
	}
}


////////////////////////////////////////////////////
// getNumRadii()
////////////////////////////////////////////////////
unsigned int BubbleData::getNumRadii() const
{
	return ( unsigned int )m_radii_list.size();
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
	std::string cmd = "select -r bubbleParticle";
	Convenience::appendNumToStdString( cmd, num );
	cmd += ";";

	MGlobal::executeCommand( Convenience::convertStdStringToMString( cmd ) );
	MGlobal::executeCommand( "doDelete" );
}