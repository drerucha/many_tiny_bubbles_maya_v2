#include "EmitterData.h"


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
void EmitterData::init( MString name )
{
	m_name = name;
}


#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>

#include "SDFGen/makelevelset3.h"

////////////////////////////////////////////////////
// generate list of positions on mesh where bubbles can emit from
// basically, fill m_source_pos_list
////////////////////////////////////////////////////
void EmitterData::setSourcePosList( const unsigned int& cubicNum,
									const unsigned int& meltingSpeed )
{
	// list of possible bubble generation locations on mesh
	std::vector<vec3> source_pos_list;

	// select Maya object by name
	MSelectionList maya_sel_list;
	MGlobal::getSelectionListByName( m_name, maya_sel_list );

	// get path to a mesh DAG node
	MDagPath dag_node_path;
	maya_sel_list.getDagPath( 0, dag_node_path );

	// if dag_node_path points to a transform node instead of a shape node, try to extend the DAG path to reach a shape node
	bool node_has_shape = true;
	if ( dag_node_path.apiType() == MFn::kTransform ) {
		MStatus stat = dag_node_path.extendToShape();
		if ( stat != MStatus::kSuccess ) {
			node_has_shape = false;
		}
	}

	// if we were able to find the shape node for the input mesh, and the node supports the kMesh function set
	if ( node_has_shape && dag_node_path.hasFn( MFn::kMesh ) ) {

		// get MFnMesh from MDagPath
		MFnMesh mesh_surface( dag_node_path );

		MFloatPointArray mesh_vertices;

		// get mesh vertices
		mesh_surface.getPoints(mesh_vertices, MSpace::kWorld);
		
		// get number of triangles in mesh and the vertex indices for each triangle
		// triangle_vertex_indices will be three times larger than triangle_num
		MIntArray triangle_num_mint_array;
		MIntArray triangle_vertex_indices_mint_array;
		mesh_surface.getTriangles( triangle_num_mint_array,
								   triangle_vertex_indices_mint_array );

		// TODO: see if we can remove this array copy operation

		// copy triangle_vertex_indices into a C++ vector
		std::vector<int> triangle_vertex_indices;
		triangle_vertex_indices.resize( triangle_vertex_indices_mint_array.length() );
		triangle_vertex_indices_mint_array.get( &triangle_vertex_indices[0] );

		// create min and max points for bounding box
		Vec3f min_bounds( std::numeric_limits<float>::max(),
						  std::numeric_limits<float>::max(),
						  std::numeric_limits<float>::max() );
		Vec3f max_bounds( -std::numeric_limits<float>::max(),
						  -std::numeric_limits<float>::max(),
						  -std::numeric_limits<float>::max());

		// use Vec3f and Vec3ui for compatibility with SDFGen level set library
		std::vector<Vec3f> vert_list;
		std::vector<Vec3ui> face_list;





		////Create face list//
		//for(int i = 0; i < triangleVertices.length(); i+=3)
		//{
		//	faceList.push_back(Vec3ui(intArray[i], intArray[i+1], intArray[i+2]));
		//}

		////Create vertices list//
		//for (int i = 0; i < points.length(); i++)
		//{
		//	Vec3f point;
		//	point[0] = points[i].x;
		//	point[1] = points[i].y;
		//	point[2] = points[i].z;
		//	vertList.push_back(point);
		//	update_minmax(point, min_box, max_box);
		//}
		////max_box is the maximum x, y, z of all the vertices 
		////min_box is the minimum x, y, z of all the vertices
		////The space between max_box and min_box could enclose all the vertices
		//	
		//double sizeX = max_box[0] - min_box[0];
		//double sizeY = max_box[1] - min_box[1];
		//double sizeZ = max_box[2] - min_box[2];
		//double maxSideLength = sizeX;
		//if(sizeY > maxSideLength)
		//	maxSideLength = sizeY;
		//if(sizeZ > maxSideLength)
		//	maxSideLength = sizeZ;
		////Here I tried to find a reasonable dx for computing Level Set Method 
		////The smaller dx will compute Level Set with finer resolution
		////The finer resolution will required more time to compute
		//double totalSize = sizeX * sizeY * sizeZ;
		//double cubicSize = totalSize/ cubicNum;
		//double dx = pow(cubicSize, (double)(1.0f/ 3.0f));

		////Add padding around the box.
		//Vec3f unit(dx, dx, dx);
		//min_box -= (float)0.001 * unit;
		//max_box += (float)0.001 * unit;
		//Vec3ui sizes = Vec3ui((max_box - min_box)/dx) + Vec3ui(2,2,2);
		//Array3f phi_grid;



		//vector<double> phiValue;
		////Use level set method to compute
		//make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);

		////phi_grid is the returning distance function
		////From returning distance function to create bubble generation source////
		//for(int i = 0 ;i < phi_grid.a.size();i++)
		//{
		//	phiValue.push_back(phi_grid.a[i]);

		//	//The bubble will generated from the outer surface of the mesh
		//	if(phi_grid.a[i] >= 0 && phi_grid.a[i] < dx/2.0f)
		//	{
		//		int z = i / phi_grid.ni / phi_grid.nj;
		//		int y = (i - z * phi_grid.ni * phi_grid.nj) / phi_grid.ni;
		//		int x = i - z * phi_grid.ni * phi_grid.nj - y * phi_grid.ni;

		//		double ran_numX = (rand()%100) / 100.0f - 0.5;
		//		double ran_numY = (rand()%100) / 100.0f - 0.5;
		//		double ran_numZ = (rand()%100) / 100.0f - 0.5;
		//		//Using random numbers to do jitter sample
		//		double gridPoxX = min_box[0] + x * dx + ran_numX * dx;
		//		double gridPosY = min_box[1] + y * dx + ran_numY * dx;
		//		double gridPosZ = min_box[2] + z * dx + ran_numZ * dx;

		//		sourcePosLocal.push_back(vec3(gridPoxX, gridPosY, gridPosZ));
		//	}
		//}

		////////////////////////////////////
		////Here is for testing melting mesh
		////If the meltingSpeed is NOT zero, than the mesh will melt
		////////////////////////////////////
		//if(meltingSpeed != 0)
		//{
		//	double* verticesValueArray = new double[phi_grid.a.size()];
		//	//After the Marching Cube Method, the points of new mesh will be stored in marchingCubePointList
		//	vector<vec3> marchingCubePointList;
		//	//After the Marching Cube Method, the face list of new mesh will be stored in triangleList
		//	vector<vector<int>> triangleList; 
		//	vector<vector<int>> temp;
		//	//The origin position is the left bottom corner of the box
		//	vec3 originPos(min_box[0], min_box[1], min_box[2]);

		//	//Compute the melting thickness, the higher value of meltingThickness means the higher speed of melting
		//	double meltingThickness = maxSideLength * meltingSpeed / 100.0f;

		//	//Convert the sign of each grid's distance function
		//	for(int i = 0 ;i < phi_grid.a.size(); i++)
		//	{
		//		verticesValueArray[i] = -phi_grid.a[i];
		//	}
		//	/////////////////////////////////////////////////////////////////////////////////
		//	//From returning distance function with marching cube method to create new mesh//
		//	/////////////////////////////////////////////////////////////////////////////////
		//	for(int z = 0 ; z < sizes[2]-1; z++)
		//	for(int y = 0 ; y < sizes[1]-1; y++)
		//	for(int x = 0 ; x < sizes[0]-1; x++)
		//	{
		//		marchingCube(meltingThickness, x, y, z,//meltingThickness
		//						sizes[0]-1, sizes[1]-1, sizes[2]-1, dx,
		//						verticesValueArray, originPos, 
		//						&marchingCubePointList,
		//						&triangleList,
		//						&temp);
		//	}
		//		
		//	////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//	//Use marchingCubePointList & triangleList to create the OBJ file and import this OBJ file as a new mesh
		//	//After importing the OBJ file, it should delete the old mesh, but it will crash MAYA with some weird reason
		//	//TODO fix the problem
		//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//	stringstream ss;
		//	string tmp;
		//	string newMeshNameStr = "meltingmesh";

		//	ss<<objectCount;
		//	ss>>tmp;
		//	newMeshNameStr += tmp;
		//	ss.clear();

		//	char* newMeshNameChar=(char *)newMeshNameStr.c_str();
		//	MString newMeshNameMStr = newMeshNameChar;
		//	string test2 = newMeshNameMStr.asChar();
		//	//MString fileName = mllPath + "/meltingmesh" +(int)objectCount + ".obj";
		//	MString fileName = mllPath + "/" + newMeshNameMStr + ".obj";
		//	string fileNameStr = fileName.asChar();
		//	createObjFile(fileNameStr, marchingCubePointList, triangleList);
		//	MString importObj = "file -import -type OBJ -ra true -namespace \"" + newMeshNameMStr + "\" -options \"mo=1\"  -pr -loadReferenceDepth \"all\" \"" + fileName + "\"";
		//	string test = importObj.asChar();
		//	MGlobal::executeCommand(importObj);
		//}





	}

	// TODO: fill m_source_pos_list directly instead of filling source_pos_list and copying

	//sourcePosLocal is used to store all the positions where generating bubbles 
	m_source_pos_list.clear();
	m_source_pos_list = source_pos_list;
}