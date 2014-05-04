//#pragma once

#ifndef _GlobalState
#define _GlobalState

#include <vector>
#include "vec.h"


class GlobalState
{

////////////////////////////////////////////////////
// methods
////////////////////////////////////////////////////

public:

	GlobalState( void );
	~GlobalState( void );

	static void setSelectedObject( std::string name );

	static void deleteSelectedObject( std::string name );

	static std::vector<std::string> getSelectedObject( void );

	static void set_maya_fluid_convert_mesh_status( bool if_mesh );

	static bool get_maya_fluid_convert_mesh_status( void );

	static void storePointList( std::vector<vec3> list );
	static std::vector<std::vector<vec3>> getStoredPointList( void );

	static void storeFaceList( std::vector<int> list );
	static std::vector<std::vector<int>> getStoredFaceList( void );

	static bool objectExists( std::string name );

	static void storeFluidTransformName( std::string name );
	static std::string getFluidTransformName( void );
};

#endif