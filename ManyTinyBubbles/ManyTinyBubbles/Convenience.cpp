#include "Convenience.h"

#include <maya/MGlobal.h>
#include <string>


////////////////////////////////////////////////////
// constructor
////////////////////////////////////////////////////
Convenience::Convenience( void )
{
}


////////////////////////////////////////////////////
// constructor
////////////////////////////////////////////////////
Convenience::~Convenience( void )
{
}


////////////////////////////////////////////////////
// std::string / MString conversions
////////////////////////////////////////////////////

std::string Convenience::convertMStringToStdString( MString mstring )
{
	std::string std_string = mstring.asChar();
	return std_string;
}

MString Convenience::convertStdStringToMString( std::string std_string )
{
	MString mstring = std_string.c_str();
	return mstring;
}


////////////////////////////////////////////////////
// print in Maya script editor
////////////////////////////////////////////////////

void Convenience::printInScriptEditor( MString mstring )
{
	MGlobal::displayInfo( mstring );
}

void Convenience::printInScriptEditor( std::string std_string )
{
	MGlobal::displayInfo( convertStdStringToMString( std_string ) );
}


////////////////////////////////////////////////////
// get parent of passed-in object
////////////////////////////////////////////////////
MString Convenience::getParent( MString child )
{
	MStringArray parent;
	MString cmd = "listRelatives -parent -path " + child;
	MGlobal::executeCommand( cmd, parent );
	return parent[0];
}


////////////////////////////////////////////////////
// append number to std::string
////////////////////////////////////////////////////
void Convenience::appendNumToStdString( std::string& str, const int& num )
{
	char num_str[21]; // large enough to hold all numbers up to 64-bits
	sprintf( num_str, "%d", num );
	str += num_str;
}


////////////////////////////////////////////////////
// attribute getters
////////////////////////////////////////////////////

int Convenience::getAttributeInt( MString object_name, MString attr_name )
{
	int attr;
	MString cmd = "getAttr " + object_name + "." + attr_name;
	MGlobal::executeCommand( cmd, attr );
	return attr;
}

double Convenience::getAttributeDouble( MString object_name, MString attr_name )
{
	double attr;
	MString cmd = "getAttr " + object_name + "." + attr_name;
	MGlobal::executeCommand( cmd, attr );
	return attr;
}

MString Convenience::getAttributeMString( MString object_name, MString attr_name )
{
	MString attr;
	MString cmd = "getAttr " + object_name + "." + attr_name;
	MGlobal::executeCommand( cmd, attr );
	return attr;
}

MIntArray Convenience::getAttributeIntArray( MString object_name, MString attr_name )
{
	MIntArray attr;
	MString cmd = "getAttr " + object_name + "." + attr_name;
	MGlobal::executeCommand( cmd, attr );
	return attr;
}

MDoubleArray Convenience::getAttributeDoubleArray( MString object_name, MString attr_name )
{
	MDoubleArray attr;
	MString cmd = "getAttr " + object_name + "." + attr_name;
	MGlobal::executeCommand( cmd, attr );
	return attr;
}

MStringArray Convenience::getAttributeMStringArray( MString object_name, MString attr_name )
{
	MStringArray attr;
	MString cmd = "getAttr " + object_name + "." + attr_name;
	MGlobal::executeCommand( cmd, attr );
	return attr;
}