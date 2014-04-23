#pragma once

#include <maya/MString.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MStringArray.h>


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

	// miscellaneous
	static MString getParent( MString child );
	static void appendNumToStdString( std::string& str, const int& num );

	// random number generation
	static int generateRandomIntInclusive( const int& min, const int& max );
	static int generateRandomIntExclusive( const int& min, const int& max );
	static double generateRandomDoubleInclusive( const double& min, const double& max );
	static double generateRandomDoubleBetweenZeroAndOneInclusive( void );
	static double generateRandomDoubleBetweenNegativeOneAndOneInclusive( void );
};