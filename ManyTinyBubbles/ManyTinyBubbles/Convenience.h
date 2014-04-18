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

	static std::string convertMStringToStdString( MString mstring );
	static MString convertStdStringToMString( std::string std_string );

	static void printInScriptEditor( MString mstring );
	static void printInScriptEditor( std::string std_string );


	////////////////////////////////////////////////////
	// attribute getters
	////////////////////////////////////////////////////

	static int getAttributeInt( MString object_name, MString attr_name );
	static double getAttributeDouble( MString object_name, MString attr_name );
	static MString getAttributeMString( MString object_name, MString attr_name );
	static MIntArray getAttributeIntArray( MString object_name, MString attr_name );
	static MDoubleArray getAttributeDoubleArray( MString object_name, MString attr_name );
	static MStringArray getAttributeMStringArray( MString object_name, MString attr_name );


	////////////////////////////////////////////////////
	// relative getters
	////////////////////////////////////////////////////

	static MString getParent( MString child );
};