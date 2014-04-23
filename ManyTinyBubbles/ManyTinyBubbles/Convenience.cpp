#include "Convenience.h"

#include <maya/MGlobal.h>

#include <string>


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////
Convenience::Convenience( void )
{
}

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
// miscellaneous
////////////////////////////////////////////////////

MString Convenience::getParent( MString child )
{
	MStringArray parent;
	MString cmd = "listRelatives -parent -path " + child;
	MGlobal::executeCommand( cmd, parent );
	return parent[0];
}

void Convenience::appendNumToStdString( std::string& str, const int& num )
{
	char num_str[21]; // large enough to hold all numbers up to 64-bits
	sprintf( num_str, "%d", num );
	str += num_str;
}


////////////////////////////////////////////////////
// random number generation
////////////////////////////////////////////////////

int Convenience::generateRandomIntInclusive( const int& min, const int& max )
{
	return rand() % ( max + 1 - min ) + min;
}

int Convenience::generateRandomIntExclusive( const int& min, const int& max )
{
	return rand() % ( max - min - 1 ) + min + 1;
}

double Convenience::generateRandomDoubleInclusive( const double& min, const double& max )
{
	//return min + static_cast<float>( rand() ) / ( static_cast<float>( RAND_MAX / ( max - min ) ) );
	return min + static_cast<double>( rand() ) / ( static_cast<double>( RAND_MAX / ( max - min ) ) );
}

double Convenience::generateRandomDoubleBetweenZeroAndOneInclusive()
{
	//return static_cast<float>( rand() ) / static_cast<float>( RAND_MAX );
	return static_cast<double>( rand() ) / static_cast<double>( RAND_MAX );
}

double Convenience::generateRandomDoubleBetweenNegativeOneAndOneInclusive()
{
	return -1.0 + static_cast<double>( rand() ) / ( static_cast<double>( RAND_MAX / ( 2.0 ) ) );
}


////////////////////////////////////////////////////
// vector indices
////////////////////////////////////////////////////

unsigned int Convenience::getIndexFromIterator( std::vector<int>::iterator it, std::vector<int> vec )
{
	return ( int )( it - vec.begin() );
}

unsigned int Convenience::getIndexFromIterator( std::vector<float>::iterator it, std::vector<float> vec )
{
	return ( int )( it - vec.begin() );
}

unsigned int Convenience::getIndexFromIterator( std::vector<double>::iterator it, std::vector<double> vec )
{
	return ( int )( it - vec.begin() );
}

unsigned int Convenience::getIndexFromIterator( std::vector<vec3>::iterator it, std::vector<vec3> vec )
{
	return ( int )( it - vec.begin() );
}

unsigned int Convenience::getIndexFromIterator( std::vector<std::vector<vec3>>::iterator it, std::vector<std::vector<vec3>> vec )
{
	return ( int )( it - vec.begin() );
}