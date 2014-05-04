#include "EmitterData.h"

#include <fstream>
#include <string>

#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>

#include "SDFGen/makelevelset3.h"
#include "Convenience.h"
#include "GlobalState.h"


////////////////////////////////////////////////////
// reference to global variable
////////////////////////////////////////////////////

extern MString mllPath;


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////

EmitterData::EmitterData()
{
}

EmitterData::~EmitterData()
{
}


////////////////////////////////////////////////////
// set m_name
////////////////////////////////////////////////////
void EmitterData::init( int emit_rate, double min_bubble_radius )
{
	m_emission_rate = emit_rate;
	m_min_bubble_radius = min_bubble_radius;
}


////////////////////////////////////////////////////
// generate list of positions on mesh where bubbles can emit from
// basically, fill m_source_pos_list
////////////////////////////////////////////////////
void EmitterData::createEmissionPositionsOnMesh( const unsigned int& voxel_num,
												 const unsigned int& melting_rate )
{
	// TODO: see if we can remove some of the data copying we're currently doing
	// TODO: fill m_source_pos_list directly instead of filling source_pos_list and copying
	// TODO: break this method up into multiple smaller methods b/c it's too long and difficult to work in

	// debug
	Convenience::printInScriptEditor( "start" );

	// list of possible bubble generation locations on mesh
	std::vector<vec3> source_pos_list;

	// get all selected mesh names
	std::vector<std::string> selectedObjectNames = GlobalState::getSelectedObject();
	for ( int objectCount = 0; objectCount< selectedObjectNames.size(); ++objectCount ) {

		// retrieve the selected mesh object
		MString strSelectObjName = ( char* )selectedObjectNames[objectCount].c_str();

		// if this object does not exist, delete the name from list
		int exist;
		MGlobal::executeCommand( "objExists " + strSelectObjName, exist );
		if ( exist == 0 ) {
			GlobalState::deleteSelectedObject( Convenience::convertMStringToStdString(strSelectObjName) );
			continue;
		}

		MStringArray emitterMeshInfoArray;

		// script to get the material name of this mesh
		MGlobal::executeCommand( "listSets -t 1 -ets -o " + strSelectObjName, emitterMeshInfoArray );
		MString emitterMeshMaterial = emitterMeshInfoArray[0];
			

		////////////////////////////////////////////////////
		// first, get Maya mesh object from mesh name
		////////////////////////////////////////////////////

		// select Maya object by name
		MSelectionList maya_sel_list;
		MGlobal::getSelectionListByName( strSelectObjName, maya_sel_list );

		// get path to a mesh DAG node
		MDagPath dag_node_path;
		maya_sel_list.getDagPath( 0, dag_node_path );

		// if dag_node_path points to a transform node instead of a shape node, try to extend the DAG path to reach a shape node
		bool node_has_shape = true;
		if ( dag_node_path.apiType() == MFn::kTransform ) {
			MStatus stat = dag_node_path.extendToShape();
			if ( stat != MStatus::kSuccess ) {
				node_has_shape = false;
				GlobalState::deleteSelectedObject( Convenience::convertMStringToStdString(strSelectObjName) );
			}
		}

		// if we were able to find the shape node for the input mesh, and the node supports the kMesh function set
		if ( node_has_shape && dag_node_path.hasFn( MFn::kMesh ) ) {

			// get MFnMesh from MDagPath
			MFnMesh mesh_surface( dag_node_path );


			////////////////////////////////////////////////////
			// next, get surface data from Maya mesh object
			////////////////////////////////////////////////////

			// get number of triangles in mesh and the vertex indices for each triangle
			// triangle_vertex_indices will be three times larger than triangle_num
			MIntArray minta_tri_num;
			MIntArray minta_tri_vertex_indices;
			mesh_surface.getTriangles( minta_tri_num,
									   minta_tri_vertex_indices );
			unsigned int triangle_vertex_num = minta_tri_vertex_indices.length();

			// copy triangle_vertex_indices into a C++ vector
			std::vector<int> triangle_vertex_indices;
			triangle_vertex_indices.resize( triangle_vertex_num );
			minta_tri_vertex_indices.get( &triangle_vertex_indices[0] );

			// create min and max points for bounding box
			// min_bounds is the minimum x, y, z of all vertices
			// max_bounds is the maximum x, y, z of all vertices
			// space enclosed between min_bounds and max_bounds enclose all vertices
			Vec3f min_bounds( std::numeric_limits<float>::max(),
							  std::numeric_limits<float>::max(),
							  std::numeric_limits<float>::max() );
			Vec3f max_bounds( -std::numeric_limits<float>::max(),
							  -std::numeric_limits<float>::max(),
							  -std::numeric_limits<float>::max() );

		
			////////////////////////////////////////////////////
			// create face list
			// use Vec3ui for compatibility with SDFGen level set library
			////////////////////////////////////////////////////

			std::vector<Vec3ui> face_list;

			// create face list
			for ( unsigned int i = 0; i < triangle_vertex_num; i += 3 ) {
				face_list.push_back( Vec3ui( triangle_vertex_indices[i],
											 triangle_vertex_indices[i + 1],
											 triangle_vertex_indices[i + 2] ) );
			}


			////////////////////////////////////////////////////
			// create vertex list
			// use Vec3f for compatibility with SDFGen level set library
			////////////////////////////////////////////////////

			std::vector<Vec3f> vert_list;

			MFloatPointArray mesh_vertices;
			mesh_surface.getPoints(mesh_vertices, MSpace::kWorld);

			for ( unsigned int i = 0; i < mesh_vertices.length(); ++i ) {
				MFloatPoint vertex = mesh_vertices[i];
				Vec3f new_vert( vertex[VX],
								vertex[VY],
								vertex[VZ] );
				vert_list.push_back( new_vert );

				// update_minmax is a method in vecL.h of SDFGen library
				// update bounding box points
				update_minmax( new_vert, min_bounds, max_bounds );
			}


			////////////////////////////////////////////////////
			// bounding box stuff
			////////////////////////////////////////////////////

			float bounding_box_size_x = max_bounds[VX] - min_bounds[VX];
			float bounding_box_size_y = max_bounds[VY] - min_bounds[VY];
			float bounding_box_size_z = max_bounds[VZ] - min_bounds[VZ];

			// get longest bounding box side
			float bounding_box_max_length = bounding_box_size_x;
			if ( bounding_box_size_y > bounding_box_max_length ) {
				bounding_box_max_length = bounding_box_size_y;
			}
			if ( bounding_box_size_z > bounding_box_max_length ) {
				bounding_box_max_length = bounding_box_size_z;
			}

			// compute reasonable dx for computing level sets
			// dx represents the level set resolution
			float bounding_box_volume = bounding_box_size_x * bounding_box_size_y * bounding_box_size_z;
			float voxel_volume = bounding_box_volume / voxel_num; // voxel_num is a method argument
			//float dx = pow( voxel_volume, ( 1.0f / 3.0f ) );
			m_level_set_dx = pow( voxel_volume, ( 1.0f / 3.0f ) );

			// add padding around bounding box
			Vec3f unit( 1.0f, 1.0f, 1.0f );
			min_bounds -= 0.001f * unit;
			max_bounds += 0.001f * unit;


			////////////////////////////////////////////////////
			// compute signed distance function
			////////////////////////////////////////////////////

			// number of sample points per x, y, z direction?
			Vec3ui sizes = Vec3ui( ( max_bounds - min_bounds ) / ( float )m_level_set_dx ) + Vec3ui( 2, 2, 2 );

			// Array3f is a data structure defined in array3.h of SDFGen library
			Array3f signed_distance_values;

			// make_level_set3 is a method in makelevelset3.cpp of SDFGen library
			make_level_set3( face_list,
							 vert_list,
							 min_bounds,
							 ( float )m_level_set_dx,
							 sizes[VX],
							 sizes[VY],
							 sizes[VZ],
							 signed_distance_values );

			// create bubble generation positions from signed distance function comuted in make_level_set3
			// signed_distance_values.a is an ArrayT
			for ( unsigned int i = 0; i < signed_distance_values.a.size(); ++i ) {

				// bubble generation position will be on outer surface of mesh or within some threshold above it
				if ( signed_distance_values.a[i] >= 0.0f && signed_distance_values.a[i] < m_level_set_dx / 2.0f ) {

					// TODO: ask what ni, nj, nk are for
					// TODO: ask about this position computation logic

					// convert linear signed distance array index into 3D indices
					int z = i / signed_distance_values.ni / signed_distance_values.nj;
					int y = ( i - z * signed_distance_values.ni * signed_distance_values.nj ) / signed_distance_values.ni;
					int x = i - z * signed_distance_values.ni * signed_distance_values.nj - y * signed_distance_values.ni;

					double gridPoxX = min_bounds[VX] + x * m_level_set_dx;
					double gridPosY = min_bounds[VY] + y * m_level_set_dx;
					double gridPosZ = min_bounds[VZ] + z * m_level_set_dx;

					// add position to source_pos_list
					source_pos_list.push_back( vec3( gridPoxX,
													 gridPosY,
													 gridPosZ ) );
				}
			}

			
			// debug
			Convenience::printInScriptEditor( "just before melting" );

			// TODO: should probably move this melting logic into its own method

			////////////////////////////////////////////////////
			// mesh melting
			////////////////////////////////////////////////////

			if ( melting_rate != 0 ) {
				double *verticesValueArray = new double[signed_distance_values.a.size()];

				// storage locations for marching cube method
				std::vector<vec3> marchingCubePointList;
				std::vector<std::vector<int>> triangleList; 
				std::vector<std::vector<int>> temp;
				
				// origin position is the left bottom corner of the bounding box
				vec3 originPos( min_bounds[0],
								min_bounds[1],
								min_bounds[2] );

				// compute melting thickness; higher value means greater melting speed
				double meltingThickness = bounding_box_max_length * melting_rate / 100.0;

				// reverse sign of each grid's distance function
				for ( unsigned int i = 0; i < signed_distance_values.a.size(); ++i ) {
					verticesValueArray[i] = -signed_distance_values.a[i];
				}


				////////////////////////////////////////////////////
				// call marching cube method
				////////////////////////////////////////////////////

				// debug
				Convenience::printInScriptEditor( "just before calling marching cube function" );

				for ( unsigned int z = 0 ; z < sizes[2] - 1; ++z ) {
					for ( unsigned int y = 0 ; y < sizes[1] - 1; ++y ) {
						for ( unsigned int x = 0 ; x < sizes[0] - 1; ++x ) {
							marchingCube( meltingThickness,
										  x, y, z,
										  sizes[0]-1, sizes[1]-1, sizes[2]-1,
										  m_level_set_dx,
										  verticesValueArray,
										  originPos,
										  &marchingCubePointList,
										  &triangleList,
										  &temp);
						}
					}
				}

				// debug
				Convenience::printInScriptEditor( "just after calling marching cube function" );


				////////////////////////////////////////////////////
				// create obj file from marchingCubePointList & triangleList
				// import this obj file as a new mesh
				// after importing, old mesh should be deleted
				////////////////////////////////////////////////////
				
				if ( bounding_box_max_length > m_min_bubble_radius ) {
					std::string newMeshNameStr;

					if ( Convenience::stringHasEnding( Convenience::convertMStringToStdString( strSelectObjName ), "_melltingmesh_Mesh" ) ) {
						std::string initialName = Convenience::convertMStringToStdString( strSelectObjName );
						newMeshNameStr = initialName.substr( 0, initialName.length() - 18 ) + "_melltingmesh";
					}
					else {
						newMeshNameStr  = Convenience::convertMStringToStdString( strSelectObjName ) + "_melltingmesh";
					}
					
					MString newMeshNameMStr = newMeshNameStr.c_str();
					
					MString fileNameMStr = mllPath + "/" + newMeshNameMStr + ".obj";
					std::string fileNameStr = fileNameMStr.asChar();

					// create obj file
					createObjFile( fileNameStr, marchingCubePointList, triangleList );

					// prepare to import this new OBJ
					MString importObj = "file -import -type OBJ -ignoreVersion -ra true -rpr \"" + newMeshNameMStr + "\" -options \"mo=1\"  -pr -loadReferenceDepth \"all\" \"" + fileNameMStr + "\"";

					// delete old mesh
					MGlobal::executeCommand( "delete " + strSelectObjName );

					// import new obj
					MGlobal::executeCommand( importObj );

					// delete the created obj file
					std::string deleteFilePath;
					if ( Convenience::stringHasEnding( Convenience::convertMStringToStdString( strSelectObjName ), "_melltingmesh_Mesh" ) ) {
						deleteFilePath = Convenience::convertMStringToStdString( mllPath + "/" + newMeshNameMStr + ".obj" );
					}
					else {
						deleteFilePath = Convenience::convertMStringToStdString( mllPath + "/" + strSelectObjName ) + "_melltingmesh.obj";
					}

					DeleteFile( deleteFilePath.c_str() );

					// delete the old stored mesh's name
					GlobalState::deleteSelectedObject( Convenience::convertMStringToStdString(strSelectObjName) );

					if ( Convenience::stringHasEnding( Convenience::convertMStringToStdString( strSelectObjName ), "_melltingmesh_Mesh" ) ) {
						GlobalState::setSelectedObject( Convenience::convertMStringToStdString( strSelectObjName ) );

						// apply shader to our new emitter mesh
						MGlobal::executeCommand( "select -r " + strSelectObjName );
						MGlobal::executeCommand( "sets -e -forceElement " + emitterMeshMaterial );
					}
					else {
						GlobalState::setSelectedObject( Convenience::convertMStringToStdString( strSelectObjName ) + "_melltingmesh_Mesh" );

						// apply shader to our new emitter mesh
						MGlobal::executeCommand( "select -r " + strSelectObjName + "_melltingmesh_Mesh" );
						MGlobal::executeCommand( "sets -e -forceElement " + emitterMeshMaterial );
					}
				}
				else {
					// debug
					Convenience::printInScriptEditor( "just before deleting old mesh" );

					// delete old mesh
					MGlobal::executeCommand( "delete " + strSelectObjName );
					GlobalState::deleteSelectedObject( Convenience::convertMStringToStdString( strSelectObjName ) );
				}

				// reclaim memory
				delete[] verticesValueArray;
			}
		}
	}

	// copy source_pos_list into the m_source_pos_list member variable
	m_source_pos_list.clear();
	m_source_pos_list = source_pos_list;

	// debug
	Convenience::printInScriptEditor( "end" );
}


////////////////////////////////////////////////////
// generate bubbles on mesh or sphere surface
// return new bubble positions, velocities, and radius group
////////////////////////////////////////////////////
void EmitterData::generateBubbles( const std::vector<double>&	bubble_radii_list,
								   std::vector<vec3>&			new_bubble_positions,
								   std::vector<vec3>&			new_bubble_velocities,
								   std::vector<unsigned int>&	new_bubble_radius_group ) const
{
	// TODO: make bubble generate rate dependent on user-defined emission rate

	int total_bubbles_per_frame = m_emission_rate;
	const vec3 INITIAL_BUBBLE_VELOCITY( 0.0, 1.0, 0.0 );

	// if emitter is a mesh
	if ( m_source_pos_list.size() != 0 ) {

		// iterate through bubble radius groups
		for ( unsigned int k = 0; k < bubble_radii_list.size(); ++k ) {

			// divide the total bubbles per frame with frame number, and round the result
			unsigned int bubbles_in_this_radii = ( int )std::floor( ( float )total_bubbles_per_frame / ( ( float )bubble_radii_list.size() - k ) + 0.5f );
			total_bubbles_per_frame -= bubbles_in_this_radii;

			if ( bubbles_in_this_radii == 0 ) {
				continue;
			}

			// if the bubbles in this radius group are more than the source postions
			if ( bubbles_in_this_radii > m_source_pos_list.size() ) {

				// iterate through every available emission position
				for ( unsigned int i = 0; i < m_source_pos_list.size(); ++i ) {

					// compute the bubbles in each position and round the number
					unsigned int bubbles_per_position = ( int )std::floor( ( float )bubbles_in_this_radii / ( float )m_source_pos_list.size() + 0.5f ) ;
					bubbles_in_this_radii -= bubbles_per_position;

					// each position has multiple bubbles, so iterate to generate bubbles
					for ( unsigned int p = 0; p < bubbles_per_position; ++p ) {
						vec3 bubble_position = m_source_pos_list[i];

						double random_num_x = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
						double random_num_y = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
						double random_num_z = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );

						// jitter the position
						bubble_position += vec3( m_level_set_dx * random_num_x,
												 m_level_set_dx * random_num_y,
												 m_level_set_dx * random_num_z );

						new_bubble_positions.push_back( bubble_position );
						new_bubble_velocities.push_back( INITIAL_BUBBLE_VELOCITY );
						new_bubble_radius_group.push_back( k );
					}
				}

			}
			// if the bubbles in this radius group are less than the source postions
			else {

				// iterate through every available emission position
				// but because the postions are more than the bubbles, bubbles are generated on period positions
				for ( unsigned int i = 0; i < bubbles_in_this_radii; ++i ) {

					//get random index from the total position list
					int position_index = Convenience::generateRandomIntInclusive( 0, ( int )m_source_pos_list.size() - 1 );
					vec3 bubble_position = m_source_pos_list[position_index];

					double random_num_x = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
					double random_num_y = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );
					double random_num_z = Convenience::generateRandomDoubleInclusive( -0.5, 0.5 );

					// jitter the position
					bubble_position += vec3( m_level_set_dx * random_num_x,
											 m_level_set_dx * random_num_y,
											 m_level_set_dx * random_num_z );

					new_bubble_positions.push_back( bubble_position );
					new_bubble_velocities.push_back( INITIAL_BUBBLE_VELOCITY );
					new_bubble_radius_group.push_back( k );
				}
			}
		}
	}
}

















////////////////////////////////////////////////////
// EmitterData::createObjFile()
////////////////////////////////////////////////////
void EmitterData::createObjFile( std::string fileName,
								 std::vector<vec3> marchingCubePointList,
								 std::vector<vector<int>> triangleList)
{
	ofstream outfile( fileName.c_str() );

	for ( unsigned int i = 0; i < marchingCubePointList.size(); ++i ) {
		outfile << "v " << marchingCubePointList[i][0] << " " << marchingCubePointList[i][1] << " " << marchingCubePointList[i][2] << std::endl;
	}

	for ( unsigned int i = 0; i < triangleList.size(); ++i ) {
		outfile << "f " << triangleList[i][0] << " " << triangleList[i][1] << " " << triangleList[i][2] << std::endl;
	}
   
	outfile.close();
}


////////////////////////////////////////////////////
// EmitterData::getSourcePosListSize()
////////////////////////////////////////////////////
int EmitterData::getSourcePosListSize()
{
	return ( int )m_source_pos_list.size();
}


////////////////////////////////////////////////////
// EmitterData::deleteTheCurrentSceneEmitterMeshs()
////////////////////////////////////////////////////
void EmitterData::deleteEmitterMeshesFromScene()
{
	// get all the selected mesh name
	std::vector<std::string> selectedObjectNames = GlobalState::getSelectedObject();

	for ( int objectCount= 0; objectCount< selectedObjectNames.size(); ++objectCount ) {

		// retrieve the selected mesh object
		MString strSelectObjName = ( char* )selectedObjectNames[objectCount].c_str();
		int exist;
		MGlobal::executeCommand( "objExists " + strSelectObjName, exist );

		if ( exist == 1 ) {
			// delete all the meshes from current scene
			MGlobal::executeCommand("delete " + strSelectObjName);
		}

		GlobalState::deleteSelectedObject( Convenience::convertMStringToStdString( strSelectObjName ) );
	}
}


////////////////////////////////////////////////////
// EmitterData::createObjFileFromStoredMeshData()
////////////////////////////////////////////////////
void EmitterData::createObjFileFromStoredMeshData()
{
	std::vector<std::vector<vec3>> pointLists = GlobalState::getStoredPointList();
	std::vector<std::vector<int>> faceLists = GlobalState::getStoredFaceList();

	for ( unsigned int i = 0; i < pointLists.size(); ++i ) {
		std::vector<vec3> pointList = pointLists[i];
		std::vector<int> faceList = faceLists[i];
		std::vector<std::vector<int>> triangleFaceList;

		// convert the linear face list index into a triangle format
		for ( unsigned int j = 0; j < faceList.size(); j += 3 ) {
			std::vector<int> newFace;
			newFace.push_back( faceList[j] );
			newFace.push_back( faceList[j + 1] );
			newFace.push_back( faceList[j + 2] );
			triangleFaceList.push_back( newFace );
		}

		std::string generalFileName = "initialMesh";
		Convenience::appendNumToStdString( generalFileName, i );
		std::string specificFileName = generalFileName + ".obj"; 
		std::string file_path_name = Convenience::convertMStringToStdString( mllPath ) + "/" + specificFileName;

		// create new mesh in the directory
		createObjFile( file_path_name, pointList, triangleFaceList );

		// store the newly created mesh name in the list of StatusData
		GlobalState::setSelectedObject(generalFileName + "_Mesh");
		std::vector<std::string> selectedObjectNames = GlobalState::getSelectedObject();
		MString importObj = "file -import -type OBJ -ignoreVersion -ra true -rpr \"" 
			+ Convenience::convertStdStringToMString( generalFileName )+ "\" -options \"mo=1\"  -pr -loadReferenceDepth \"all\" \"" 
			+ Convenience::convertStdStringToMString( file_path_name ) + "\"";

		// import new obj
		MGlobal::executeCommand( importObj );

		// delete the created OBJ file
		DeleteFile( file_path_name.c_str() );
	}
}


////////////////////////////////////////////////////
// EmitterData::vertexInterp() - used in marching cube method
////////////////////////////////////////////////////
vec3 EmitterData::vertexInterp( double	&isolevel,
								vec3	p1,
								vec3	p2,
								double	&valueP1,
								double	&valueP2)
{
	vec3 interPt;
	double mu;

	if ( abs( isolevel - valueP1 ) < 0.00001 ) {
		return p1;
	}
	if ( abs( isolevel - valueP2 ) < 0.00001 ) {
		return p2;
	}
	if ( abs( valueP1 - valueP2 ) < 0.00001 ) {
		return p1;
	}

	mu = ( isolevel - valueP1 ) / ( valueP2 - valueP1 );
	interPt[0] = p1[0] + mu * (p2[0] - p1[0]);
	interPt[1] = p1[1] + mu * (p2[1] - p1[1]);
	interPt[2] = p1[2] + mu * (p2[2] - p1[2]);

	return interPt;
}


////////////////////////////////////////////////////
// EmitterData::marchingCube()
////////////////////////////////////////////////////
int EmitterData::marchingCube( double							&isolevel,
							   const unsigned int				&indexX,
							   const unsigned int				&indexY,
							   const unsigned int				&indexZ,
							   const unsigned int				&meshResX,
							   const unsigned int				&meshResY,
							   const unsigned int				&meshResZ,
							   double							&meshSize,
							   double							*verticesValueArray,
							   vec3								originPos,
							   std::vector<vec3>				*marchingCubePointList,
							   std::vector<std::vector<int>>	*triangleList,
							   std::vector<std::vector<int>>	*temp)
{
	int gridValueIndex0 = indexZ * (meshResX + 1) * (meshResY + 1) + indexY * (meshResX + 1) + indexX;
	int gridValueIndex1 = indexZ * (meshResX + 1) * (meshResY + 1) + indexY * (meshResX + 1) + indexX + 1;
	int gridValueIndex2 = indexZ * (meshResX + 1) * (meshResY + 1) + (indexY + 1) * (meshResX + 1) + indexX + 1;
	int gridValueIndex3 = indexZ * (meshResX + 1) * (meshResY + 1) + (indexY + 1) * (meshResX + 1) + indexX;
	int gridValueIndex4 = (indexZ + 1) * (meshResX + 1) * (meshResY + 1) + indexY * (meshResX + 1) + indexX;
	int gridValueIndex5 = (indexZ + 1) * (meshResX + 1) * (meshResY + 1) + indexY * (meshResX + 1) + indexX + 1;
	int gridValueIndex6 = (indexZ + 1) * (meshResX + 1) * (meshResY + 1) + (indexY + 1) * (meshResX + 1) + indexX + 1;
	int gridValueIndex7 = (indexZ + 1) * (meshResX + 1) * (meshResY + 1) + (indexY + 1) * (meshResX + 1) + indexX;

	double gridValue0 = verticesValueArray[gridValueIndex0];
	double gridValue1 = verticesValueArray[gridValueIndex1];
	double gridValue2 = verticesValueArray[gridValueIndex2];
	double gridValue3 = verticesValueArray[gridValueIndex3];
	double gridValue4 = verticesValueArray[gridValueIndex4];
	double gridValue5 = verticesValueArray[gridValueIndex5];
	double gridValue6 = verticesValueArray[gridValueIndex6];
	double gridValue7 = verticesValueArray[gridValueIndex7];

	vec3 gridPos0(indexX       * meshSize + originPos[0], indexY       * meshSize + originPos[1], indexZ       * meshSize + originPos[2]);
	vec3 gridPos1((indexX + 1) * meshSize + originPos[0], indexY       * meshSize + originPos[1], indexZ       * meshSize + originPos[2]);
	vec3 gridPos2((indexX + 1) * meshSize + originPos[0], (indexY + 1) * meshSize + originPos[1], indexZ       * meshSize + originPos[2]);
	vec3 gridPos3(indexX       * meshSize + originPos[0], (indexY + 1) * meshSize + originPos[1], indexZ       * meshSize + originPos[2]);
	vec3 gridPos4(indexX       * meshSize + originPos[0], indexY       * meshSize + originPos[1], (indexZ + 1) * meshSize + originPos[2]);
	vec3 gridPos5((indexX + 1) * meshSize + originPos[0], indexY       * meshSize + originPos[1], (indexZ + 1) * meshSize + originPos[2]);
	vec3 gridPos6((indexX + 1) * meshSize + originPos[0], (indexY + 1) * meshSize + originPos[1], (indexZ + 1) * meshSize + originPos[2]);
	vec3 gridPos7(indexX       * meshSize + originPos[0], (indexY + 1) * meshSize + originPos[1], (indexZ + 1) * meshSize + originPos[2]);


	////////////////////////////////////////////////////
	// create list
	////////////////////////////////////////////////////

	int edgeTable[256]={
		0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
		0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
		0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
		0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
		0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
		0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
		0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
		0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
		0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
		0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
		0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
		0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
		0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
		0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
		0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
		0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
		0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
		0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
		0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
		0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
		0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
		0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
		0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
		0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
		0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
		0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
		0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
		0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
		0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
		0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
		0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
		0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };


	int triTable[256][16] =
		{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
		{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
		{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
		{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
		{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
		{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
		{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
		{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
		{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
		{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
		{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
		{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
		{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
		{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
		{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
		{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
		{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
		{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
		{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
		{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
		{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
		{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
		{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
		{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
		{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
		{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
		{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
		{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
		{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
		{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
		{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
		{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
		{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
		{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
		{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
		{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
		{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
		{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
		{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
		{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
		{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
		{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
		{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
		{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
		{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
		{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
		{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
		{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
		{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
		{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
		{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
		{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
		{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
		{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
		{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
		{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
		{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
		{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
		{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
		{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
		{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
		{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
		{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
		{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
		{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
		{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
		{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
		{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
		{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
		{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
		{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
		{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
		{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
		{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
		{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
		{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
		{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
		{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
		{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
		{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
		{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
		{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
		{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
		{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
		{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
		{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
		{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
		{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
		{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
		{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
		{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
		{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
		{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
		{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
		{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
		{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
		{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
		{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
		{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
		{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
		{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
		{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
		{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
		{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
		{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
		{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
		{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
		{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
		{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
		{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
		{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
		{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
		{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
		{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
		{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
		{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
		{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
		{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
		{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
		{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
		{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
		{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
		{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
		{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
		{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
		{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
		{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
		{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
		{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
		{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
		{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
		{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
		{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
		{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
		{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
		{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
		{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
		{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
		{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
		{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
		{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
		{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
		{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
		{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
		{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
		{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
		{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
		{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
		{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
		{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
		{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
		{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
		{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
		{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
		{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
		{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
		{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
		{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
		{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
		{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
		{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
		{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
		{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
		{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
		{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
		{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
		{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
		{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
		{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
		{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
		{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
		{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
		{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
		{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
		{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
		{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
		{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
		{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
		{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
		{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
		{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
	

	// determine the index into the edge table which tells us which vertices are inside the surface
	
	int cubeindex = 0;
	if (gridValue0 > isolevel) cubeindex |= 1;
	if (gridValue1 > isolevel) cubeindex |= 2;
	if (gridValue2 > isolevel) cubeindex |= 4;
	if (gridValue3 > isolevel) cubeindex |= 8;
	if (gridValue4 > isolevel) cubeindex |= 16;
	if (gridValue5 > isolevel) cubeindex |= 32;
	if (gridValue6 > isolevel) cubeindex |= 64;
	if (gridValue7 > isolevel) cubeindex |= 128;

	// Cube is entirely in/out of the surface //
	std::vector<int> pointList_eachgrid;
	if (edgeTable[cubeindex] == 0)
	{
		temp->push_back(pointList_eachgrid);
		return(0);
	}


	////////////////////////////////////////////////////
	// find the vertices where the surface intersects the cube
	////////////////////////////////////////////////////

	int edge[12];
	if (edgeTable[cubeindex] & 1)//EDGE 0
	{
		if(indexY == 0 && indexZ == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos0, gridPos1, gridValue0, gridValue1);
			//vec3 vertex = (gridPos0 + gridPos1) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[0] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[0]);

		}
		else if(indexY == 0)
		{
			int vertexIndex = (indexZ-1) * meshResX  * meshResY + indexX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[0] = temp2[4];
			//edge[0] = temp[vertexIndex][4];
			pointList_eachgrid.push_back(edge[0]);
		}
		else if(indexZ == 0)
		{
			int vertexIndex = (indexY-1) * meshResX + indexX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[0] = temp2[2];
			//edge[0] = temp[vertexIndex][2];
			pointList_eachgrid.push_back(edge[0]);
		}
		else
		{
			int vertexIndex = (indexZ-1) * meshResX  * meshResY + (indexY-1) * meshResX + indexX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[0] = temp2[6];
			//edge[0] = temp[vertexIndex][6];
			pointList_eachgrid.push_back(edge[0]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 2)//EDGE 1
	{
		if(indexZ == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos1, gridPos2, gridValue1, gridValue2);
			//vec3 vertex = (gridPos1 + gridPos2) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[1] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[1]);
		}
		else
		{
			int vertexIndex = (indexZ-1) * meshResX  * meshResY + indexY * meshResX + indexX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[1] = temp2[5];
			//edge[1] = temp[vertexIndex][5];
			pointList_eachgrid.push_back(edge[1]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 4)//EDGE 2
	{
		if(indexZ == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos2, gridPos3, gridValue2, gridValue3);
			//vec3 vertex = (gridPos2 + gridPos3) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[2] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[2]);
		}
		else
		{
			int vertexIndex = (indexZ-1) * meshResX  * meshResY + indexY * meshResX + indexX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[2] = temp2[6];
			//edge[2] = temp[vertexIndex][6];
			pointList_eachgrid.push_back(edge[2]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 8)//EDGE 3
	{
		if(indexX == 0 && indexZ == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos3, gridPos0, gridValue3, gridValue0);
			//vec3 vertex = (gridPos3 + gridPos0) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[3] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[3]);
		}
		else if(indexX == 0)
		{
			int vertexIndex = (indexZ-1) * meshResX  * meshResY + indexY * meshResX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[3] = temp2[7];
			//edge[3] = temp[vertexIndex][7];
			pointList_eachgrid.push_back(edge[3]);
		}
		else if(indexZ == 0)
		{
			int vertexIndex = indexY * meshResX + indexX-1;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[3] = temp2[1];
			//edge[3] = temp[vertexIndex][1];
			pointList_eachgrid.push_back(edge[3]);
		}
		else
		{
			int vertexIndex = (indexZ-1) * meshResX  * meshResY + indexY * meshResX + indexX-1;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[3] = temp2[5];
			//edge[3] = temp[vertexIndex][5];
			pointList_eachgrid.push_back(edge[3]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 16)//EDGE 4
	{
		if(indexY == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos4, gridPos5, gridValue4, gridValue5);
			//vec3 vertex = (gridPos4 + gridPos5) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[4] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[4]);
		}
		else
		{
			int vertexIndex = indexZ * meshResX  * meshResY + (indexY-1) * meshResX + indexX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[4] = temp2[6];
			//edge[4] = temp[vertexIndex][6];
			pointList_eachgrid.push_back(edge[4]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 32)//EDGE 5
	{
		vec3 vertex = vertexInterp(isolevel, gridPos5, gridPos6, gridValue5, gridValue6);
		//vec3 vertex = (gridPos5 + gridPos6) / 2.0f;
		marchingCubePointList->push_back(vertex);
		edge[5] = ( int )marchingCubePointList->size();
		pointList_eachgrid.push_back(edge[5]);
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 64)//EDGE 6
	{
		vec3 vertex = vertexInterp(isolevel, gridPos6, gridPos7, gridValue6, gridValue7);
		//vec3 vertex = (gridPos6 + gridPos7) / 2.0f;
		marchingCubePointList->push_back(vertex);
		edge[6] = ( int )marchingCubePointList->size();
		pointList_eachgrid.push_back(edge[6]);
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 128)//EDGE 7
	{
		if(indexX == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos7, gridPos4, gridValue7, gridValue4);
			//vec3 vertex = (gridPos7 + gridPos4) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[7] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[7]);
		}
		else
		{
			int vertexIndex = indexZ * meshResX  * meshResY + indexY * meshResX + indexX-1;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[7] = temp2[5];
			//edge[7] = temp[vertexIndex][5];
			pointList_eachgrid.push_back(edge[7]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 256)//EDGE 8
	{
		if(indexX == 0 && indexY == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos0, gridPos4, gridValue0, gridValue4);
			//vec3 vertex = (gridPos0 + gridPos4) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[8] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[8]);
		}
		else if(indexX == 0)
		{
			int vertexIndex = indexZ * meshResX  * meshResY + (indexY-1) * meshResX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[8] = temp2[11];
			//edge[8] = temp[vertexIndex][11];
			pointList_eachgrid.push_back(edge[8]);
		}
		else if(indexY == 0)
		{
			int vertexIndex = indexZ * meshResX  * meshResY + indexX-1;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[8] = temp2[9];
			//edge[8] = temp[vertexIndex][9];
			pointList_eachgrid.push_back(edge[8]);
		}
		else
		{
			int vertexIndex = indexZ * meshResX  * meshResY + (indexY-1) * meshResX + indexX-1;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[8] = temp2[10];
			//edge[8] = temp[vertexIndex][10];
			pointList_eachgrid.push_back(edge[8]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 512)//EDGE 9
	{
		if(indexY == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos1, gridPos5, gridValue1, gridValue5);
			//vec3 vertex = (gridPos1 + gridPos5) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[9] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[9]);
		}
		else
		{
			int vertexIndex = indexZ * meshResX  * meshResY + (indexY-1) * meshResX + indexX;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[9] = temp2[10];
			//edge[9] = temp[vertexIndex][10];
			pointList_eachgrid.push_back(edge[9]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 1024)//EDGE 10
	{
		vec3 vertex = vertexInterp(isolevel, gridPos2, gridPos6, gridValue2, gridValue6);
		//vec3 vertex = (gridPos2 + gridPos6) / 2.0f;
		marchingCubePointList->push_back(vertex);
		edge[10] = ( int )marchingCubePointList->size();
		pointList_eachgrid.push_back(edge[10]);
	}
	else
		pointList_eachgrid.push_back(-1);

	if (edgeTable[cubeindex] & 2048)//EDGE 11
	{
		if(indexX == 0)
		{
			vec3 vertex = vertexInterp(isolevel, gridPos3, gridPos7, gridValue3, gridValue7);
			//vec3 vertex = (gridPos3 + gridPos7) / 2.0f;
			marchingCubePointList->push_back(vertex);
			edge[11] = ( int )marchingCubePointList->size();
			pointList_eachgrid.push_back(edge[11]);
		}
		else
		{
			int vertexIndex = indexZ * meshResX  * meshResY + indexY * meshResX + indexX-1;
			vector<int> temp2 = temp->at(vertexIndex);
			edge[11] = temp2[10];
			//edge[11] = temp[vertexIndex][10];
			pointList_eachgrid.push_back(edge[11]);
		}
	}
	else
		pointList_eachgrid.push_back(-1);

	temp->push_back(pointList_eachgrid);

	
	////////////////////////////////////////////////////
	// create the triangle
	////////////////////////////////////////////////////

	for (int i = 0; triTable[cubeindex][i] != -1; i+=3) 
	{
		vector<int> triangle;
		triangle.push_back(edge[triTable[cubeindex][i  ]]);
		triangle.push_back(edge[triTable[cubeindex][i+1]]);
		triangle.push_back(edge[triTable[cubeindex][i+2]]);

		triangleList->push_back(triangle);
	}


	return 0;
}