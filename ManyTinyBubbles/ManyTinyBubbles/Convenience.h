#pragma once

#include <maya/MString.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MStringArray.h>
//#include <maya/MDagPath.h>

#include <vector>
#include "vec.h"


class Convenience
{
public:
	Convenience( void );
	~Convenience (void );

	// std::string / MString conversion
	static std::string convertMStringToStdString( MString mstring );
	static MString convertStdStringToMString( std::string std_string );

	// print in Maya's script editor
	static void printInScriptEditor( MString mstring );
	static void printInScriptEditor( std::string std_string );

	// attribute getters
	static int getAttributeInt( MString object_name, MString attr_name );
	static double getAttributeDouble( MString object_name, MString attr_name );
	static MString getAttributeMString( MString object_name, MString attr_name );
	static MIntArray getAttributeIntArray( MString object_name, MString attr_name );
	static MDoubleArray getAttributeDoubleArray( MString object_name, MString attr_name );
	static MStringArray getAttributeMStringArray( MString object_name, MString attr_name );

	// attribute setters
	//static void setAttributeInt( const MString&	object_name,
	//							 const MString&	attr_name,
	//							 const int&		val );

	// miscellaneous
	static MString getParent( MString child );
	static void setParticleRenderTypeToSphere( const MString& particle_name );
	//static bool meshNameDoesCorrespondToMayaMeshObject( const MString& mesh_name );
	//static MDagPath getDagPathToMeshNodeFromName( const MString& mesh_name );

	// std::string number concatenation
	static void appendNumToStdString( std::string& str, const int& num );
	static void appendNumToStdString( std::string& str, const unsigned int& num );
	static void appendNumToStdString( std::string& str, const float& num );
	static void appendNumToStdString( std::string& str, const double& num );
	static void appendVec3ToStdString( std::string& str, const vec3& vec );

	// random number generation
	static int generateRandomIntInclusive( const int& min, const int& max );
	static int generateRandomIntExclusive( const int& min, const int& max );
	static double generateRandomDoubleInclusive( const double& min, const double& max );
	static double generateRandomDoubleBetweenZeroAndOneInclusive( void );
	static double generateRandomDoubleBetweenNegativeOneAndOneInclusive( void );

	// vector indices
	static unsigned int getIndexFromIterator( std::vector<int>::iterator it, std::vector<int> vec );
	static unsigned int getIndexFromIterator( std::vector<float>::iterator it, std::vector<float> vec );
	static unsigned int getIndexFromIterator( std::vector<double>::iterator it, std::vector<double> vec );
	static unsigned int getIndexFromIterator( std::vector<vec3>::iterator it, std::vector<vec3> vec );
	static unsigned int getIndexFromIterator( std::vector<std::vector<vec3>>::iterator it, std::vector<std::vector<vec3>> vec );
};