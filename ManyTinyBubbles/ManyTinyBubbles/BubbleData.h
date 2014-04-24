#pragma once

#include "vec.h"
#include <vector>


// TODO: make NUM_RADII a user-definable node attribute

// constants
const int NUM_RADII = 10;

class BubbleData
{

////////////////////////////////////////////////////
// methods
////////////////////////////////////////////////////

public:

	BubbleData( void );
	~BubbleData( void );

	void init( double scattering_frequency,
			   double scattering_coefficient,
			   double breakup_frequency,
			   double size_min,
			   double size_max );

	void deleteAllParticlesInMaya();
	void reset();

	// int getters
	unsigned int getNumRadii( void ) const;
	
	// double getters
	double getRadiusAtIndex( const unsigned int& index ) const;
	double getScatteringFrequency( void ) const;
	double getScatteringCoefficient( void ) const;
	double getBreakupFrequency( void ) const;

	// std::vector getters
	std::vector<std::vector<vec3>> getPosList( void ) const;
	std::vector<double> getRadiiList( void ) const;
	//std::vector<vec3> getPosListForRadiusGroupAtIndex( const unsigned int& i ) const;

	void removeBubbleAtIndex( const unsigned int& i,
							  const unsigned int& j );
	void addBubblePosToRadiusGroupAtIndex( const vec3& pos,
										   const unsigned int& i );

	void createMayaParticlesWithName( const std::string& particle_name ) const;


private:

	void setRadii( const double&	radius_min,
				   const double&	radius_max );

	int checkIfParticleExists( const unsigned int& num ) const;
	void deleteParticle( const unsigned int& num ) const;


////////////////////////////////////////////////////
// members
////////////////////////////////////////////////////

private:
	
	// TODO: find out what m_vel_list is used for

	std::vector<std::vector<vec3>> m_pos_list;
	std::vector<std::vector<vec3>> m_vel_list;
	std::vector<double> m_radii_list;

	double m_scattering_frequency;
	double m_scattering_coefficient;
	double m_breakup_frequency;
	double m_size_min;
	double m_size_max;

};