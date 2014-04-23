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

	// getters
	unsigned int getNumRadii( void ) const;
	std::vector<std::vector<vec3>> getPosList( void ) const;
	double getRadiusAtIndex( unsigned int& index ) const;
	double getScatteringFrequency( void ) const;
	double getScatteringCoefficient( void ) const;

	void removeBubbleAtIndex( const unsigned int& i,
							  const unsigned int& j );

private:

	void setRadii( const double&	radius_min,
				   const double&	radius_max );

	int checkIfParticleExists( const unsigned int& num ) const;
	void deleteParticle( const unsigned int& num ) const;


////////////////////////////////////////////////////
// members
////////////////////////////////////////////////////

private:

	std::vector<std::vector<vec3>> m_pos_list;
	std::vector<std::vector<vec3>> m_vel_list;
	std::vector<double> m_radii_list;

	double m_scattering_frequency;
	double m_scattering_coefficient;
	double m_breakup_frequency;
	double m_size_min;
	double m_size_max;

};