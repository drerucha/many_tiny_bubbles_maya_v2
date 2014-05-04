#include "GlobalState.h"

#include <string>


////////////////////////////////////////////////////
// member variables
////////////////////////////////////////////////////

static std::vector<std::string> m_select_object_array;
static bool m_fluid_mesh_status;
static std::vector<std::vector<vec3>> m_initial_mesh_point_list;
static std::vector<std::vector<int>> m_initial_mesh_triangle_list;
static std::string m_fluid_transform_name;


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

GlobalState::GlobalState()
{
}

GlobalState::~GlobalState()
{
}


///////////////////////////////////////////////////////////////////
// when selecting meshes to be the bubble source, store the mesh's name in m_select_object_array
///////////////////////////////////////////////////////////////////
void GlobalState::setSelectedObject( std::string name )
{
	// if the name has been stored in the selectObjectArray, then skip it
	for ( unsigned int i = 0; i < m_select_object_array.size(); ++i ) {
		std::string nameOnList = m_select_object_array[i];
		if ( nameOnList == name ) {
			return;
		}
	}
	m_select_object_array.push_back( name );
}


///////////////////////////////////////////////////////////////////
// when the mesh is no longer exist, delete it from list
///////////////////////////////////////////////////////////////////
void GlobalState::deleteSelectedObject( std::string name )
{
	int i = 0;
	for ( std::vector<std::string>::iterator iter = m_select_object_array.begin(); iter != m_select_object_array.end(); ++i, ++iter ) {
		std::string nameOnList = m_select_object_array[i];
		if ( nameOnList == name ) {
			m_select_object_array.erase( iter );
			break;
		}
	}
}


///////////////////////////////////////////////////////////////////
// return all the stored mesh names
///////////////////////////////////////////////////////////////////
std::vector<std::string> GlobalState::getSelectedObject()
{
	return m_select_object_array;
}


///////////////////////////////////////////////////////////////////
// get / set the simulation processing flag
///////////////////////////////////////////////////////////////////

bool GlobalState::get_maya_fluid_convert_mesh_status()
{
	return m_fluid_mesh_status;
}

void GlobalState::set_maya_fluid_convert_mesh_status( bool if_mesh )
{
	m_fluid_mesh_status = if_mesh;
}


///////////////////////////////////////////////////////////////////
// set point list of passed-in mesh emitter
///////////////////////////////////////////////////////////////////
void GlobalState::storePointList( std::vector<vec3> list )
{
	m_initial_mesh_point_list.push_back( list );
}


///////////////////////////////////////////////////////////////////
// get point lists for all mesh emitters
///////////////////////////////////////////////////////////////////
std::vector<std::vector<vec3>> GlobalState::getStoredPointList()
{
	return m_initial_mesh_point_list;
}


///////////////////////////////////////////////////////////////////
// set face list of passed-in mesh emitter
///////////////////////////////////////////////////////////////////
void GlobalState::storeFaceList( std::vector<int> list )
{
	m_initial_mesh_triangle_list.push_back( list );
}


///////////////////////////////////////////////////////////////////
// get face lists for all mesh emitters
///////////////////////////////////////////////////////////////////
std::vector<std::vector<int>> GlobalState::getStoredFaceList()
{
	return m_initial_mesh_triangle_list;
}


///////////////////////////////////////////////////////////////////
// determine whether the object exists or not
///////////////////////////////////////////////////////////////////
bool GlobalState::objectExists( std::string name )
{
	// if name is stored in m_select_object_array, then return true
	for ( unsigned int i = 0; i < m_select_object_array.size(); ++i ) {
		std::string nameOnList = m_select_object_array[i];
		if ( nameOnList == name ) {
			return true;
		}
	}
	
	return false;
}


///////////////////////////////////////////////////////////////////
// get / set Maya fluid transform name
///////////////////////////////////////////////////////////////////

void GlobalState::storeFluidTransformName( std::string name)
{
	m_fluid_transform_name = name;
}

std::string GlobalState::getFluidTransformName()
{
	return m_fluid_transform_name;
}