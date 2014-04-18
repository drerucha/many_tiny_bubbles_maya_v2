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


////////////////////////////////////////////////////
// relative getters
////////////////////////////////////////////////////

MString Convenience::getParent( MString child )
{
	MString parent;
	MString cmd = "listRelatives -parent -path " + child;
	MGlobal::executeCommand( cmd, parent );
	return parent;
}