//#pragma once

#ifndef _BubbleData
#define _BubbleData

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

	// getters
	unsigned int					getNumRadii( void ) const;
	double							getRadiusAtIndex( const unsigned int& index ) const;
	double							getScatteringFrequency( void ) const;
	double							getScatteringCoefficient( void ) const;
	double							getBreakupFrequency( void ) const;
	std::vector<std::vector<vec3>>	getPosList( void ) const;
	std::vector<double>				getRadiiList( void ) const;
	vec3							getVelocityAtIndex( const unsigned int& i, const unsigned int& j ) const;
	//std::vector<vec3>				getPosListForRadiusGroupAtIndex( const unsigned int& i ) const;
	unsigned int					getTotalBubbleNumber() const;

	void addBubblePosToRadiusGroupAtIndex( const vec3& pos, const unsigned int& i );

	void setVelocityAtIndex( const vec3& new_vel,
							 const unsigned int& i,
							 const unsigned int& j );

	//void updatePosAtIndex( const unsigned int& i, const unsigned int& j );

	void updateBubblePositionsAtIndex( const unsigned int&	i,
									   const unsigned int&	j,
									   const float&			time_step );

	void init( double scattering_frequency,
			   double scattering_coefficient,
			   double breakup_frequency,
			   double size_min,
			   double size_max ); // Danny was here
	void reset();  // Danny was here
	void updateBubblePositions( const float& time_step ); // Danny was here
	void createMayaParticlesWithName( const std::string& particle_name ) const; // Danny was here
	void setRadiiForMayaParticlesWithName( const std::string& particle_name ) const; // Danny was here
	void deleteAllParticlesInMaya(); // Danny was here
	void addBubblePosToRadiusGroupAtIndex( std::vector<vec3> pos_list, std::vector<unsigned int> radius_group_index_list ); // Danny was here
	void addBubbleVelToRadiusGroupAtIndex( std::vector<vec3> vel_list, std::vector<unsigned int> radius_group_index_list ); // Danny was here
	unsigned int getNumRadiusGroups( void ) const; // Danny was here
	unsigned int getNumBubblesInListWithIndex( const unsigned int& i ) const; // Danny was here
	void removeBubbleAtIndex( const unsigned int& i, const unsigned int& j ); // Danny was here
	vec3 getPosAtIndex( const unsigned int& i, const unsigned int& j ) const; // Danny was here
	vec3 getVelAtIndex( const unsigned int& i, const unsigned int& j ) const; // Danny was here
	void addBubble( const unsigned int&	radius_group_index,
					vec3				pos,
					vec3				vel ); // Danny was here

private:

	void setRadii( const double& radius_min, const double& radius_max );

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

#endif