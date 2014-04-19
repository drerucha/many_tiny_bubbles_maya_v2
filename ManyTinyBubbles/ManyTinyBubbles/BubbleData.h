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

	void init( float scattering_frequency,
			   float scattering_coefficient,
			   float breakup_frequency,
			   float size_min,
			   float size_max );

	void deleteAllParticlesInMaya();
	void reset();

	// getters
	unsigned int getNumRadii( void ) const;
	std::vector<std::vector<vec3>> getPosList( void ) const;

private:

	void setRadii( const float&	radius_min,
				   const float&	radius_max );

	int checkIfParticleExists( const unsigned int& num ) const;
	void deleteParticle( const unsigned int& num ) const;


////////////////////////////////////////////////////
// members
////////////////////////////////////////////////////

private:

	std::vector<std::vector<vec3>> m_pos_list;
	std::vector<std::vector<vec3>> m_vel_list;
	std::vector<float> m_radii_list;

	float m_scattering_frequency;
	float m_scattering_coefficient;
	float m_breakup_frequency;
	float m_size_min;
	float m_size_max;

};