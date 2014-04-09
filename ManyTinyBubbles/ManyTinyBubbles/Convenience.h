#pragma once

#include <maya/MString.h>


class Convenience
{
public:
	Convenience( void );
	~Convenience (void );

	static std::string convertMStringToStdString( MString mstring );
	static MString convertStdStringToMString( std::string std_string );

	static void printInScriptEditor( MString mstring );
	static void printInScriptEditor( std::string std_string );
};

