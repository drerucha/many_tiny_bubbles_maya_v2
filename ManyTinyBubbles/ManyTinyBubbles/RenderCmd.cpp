#include "RenderCmd.h"

#include <sstream>
#include <direct.h>

#include <maya/MSyntax.h>
#include <maya/MGlobal.h>
#include <maya/MArgDatabase.h>

#include "Convenience.h"
#include "GlobalState.h"


////////////////////////////////////////////////////
// reference to global variable
////////////////////////////////////////////////////

extern MString mllPath;


////////////////////////////////////////////////////
// constants
////////////////////////////////////////////////////

// short flag can't be longer than 3 char
// long flag can't be less than 4 char
const char *start_frame_flag_short	= "-sf",	*start_frame_flag_long	= "-start_frame";
const char *end_frame_flag_short	= "-ef",	*end_frame_flag_long	= "-end_frame";


////////////////////////////////////////////////////
// constructor / destructor
////////////////////////////////////////////////////
RenderCmd::RenderCmd() : MPxCommand()
{
}

RenderCmd::~RenderCmd()
{
}


////////////////////////////////////////////////////
// syntax()
////////////////////////////////////////////////////
MSyntax RenderCmd::syntax()
{
    MSyntax syntax;

	syntax.addFlag( start_frame_flag_short,
					start_frame_flag_long,
					MSyntax::kLong );
	syntax.addFlag( end_frame_flag_short,
					end_frame_flag_long,
					MSyntax::kLong );

	return syntax;
}


////////////////////////////////////////////////////
// doIt()
////////////////////////////////////////////////////
MStatus RenderCmd::doIt( const MArgList& args )
{
	// TODO: remove "polySurface1" hard-coded name
	// TODO: make a convenience command to convert Maya fluid to polygons given a name


	////////////////////////////////////////////////////
	// convert Maya fluid to polygon mesh
	////////////////////////////////////////////////////

	int fluid_exists;
	MGlobal::executeCommand( "objExists polySurface1", fluid_exists );

	if ( fluid_exists == 0 ) {

		// TODO: uncomment the following line after creating StatusData class

		// convert Maya fluid to polygon mesh
		std::string fluid_transform_name = GlobalState::getFluidTransformName();
		MGlobal::executeCommand( "select -r " + Convenience::convertStdStringToMString( fluid_transform_name ) );
		MGlobal::executeCommand( "fluidToPoly" );
	}


	////////////////////////////////////////////////////
	// attach water shader
	////////////////////////////////////////////////////

	int shader_exists;
	MGlobal::executeCommand( "objExists seasky_water_waterBlinnSG", shader_exists );

	if ( shader_exists == 0 ) {
		// import shader
		MString fileNameMStr = mllPath + "/seasky.mb";
		MString importShader = "file -import -type \"mayaBinary\" -ignoreVersion -ra true -rpr \"seasky\" -options \"v=0;p=17\"  -pr -loadReferenceDepth \"all\" \"" + fileNameMStr + "\"";
		MGlobal::executeCommand( importShader );
		// delete the unecessary floor and surround sphere
		MGlobal::executeCommand( "delete seasky_seaPlane1" );
		MGlobal::executeCommand( "delete seasky_skySphere3" );
	}

	// TODO: fix this part

	// get the material of an object
	MStringArray waterMeshInfoArray;
	MGlobal::executeCommand( "listSets -t 1 -ets -o polySurface1", waterMeshInfoArray );
	MString waterMeshMaterial = waterMeshInfoArray[0];
	if ( waterMeshMaterial != "seasky_water_waterBlinnSG" ) {
		// apply the shader to our water polygon mesh
		MGlobal::executeCommand( "select -r polySurface1" );
		MGlobal::executeCommand( "sets -e -forceElement seasky_water_waterBlinnSG" );
	}


	////////////////////////////////////////////////////
	// get command parameters
	////////////////////////////////////////////////////

	MString name;
	MString id;

	MGlobal::displayInfo( "Render Sequential Images" );

	unsigned int start_frame = 1;
	unsigned int end_frame = 1;

	MArgDatabase argData( syntax(), args );

	if ( argData.isFlagSet( start_frame_flag_short ) ) {
		argData.getFlagArgument( start_frame_flag_short, 0, start_frame );
	}
	if ( argData.isFlagSet( end_frame_flag_short ) ) {
		argData.getFlagArgument( end_frame_flag_short, 0, end_frame );
	}

	// TODO: add warning message

	if ( start_frame < 1 ) {
		start_frame = 1;
	}
	if ( end_frame < start_frame ) {
		end_frame = start_frame;
	}


	////////////////////////////////////////////////////
	// render simulation
	////////////////////////////////////////////////////

	// set render type to mental ray
	std::string renderType = ( std::string )"setAttr " + '\"' + ( std::string )"defaultRenderGlobals.currentRenderer" + '\"' + ( std::string )" -type " + '\"' + ( std::string )"string" + '\"' +" "+ '\"'+ ( std::string )"mentalRay"+ '\"';
	char *renderTypeChar = ( char* )renderType.c_str();
	MString strRenderType = renderTypeChar;
	MGlobal::executeCommand( strRenderType );

	MGlobal::executeCommand( "setAttr defaultRenderGlobals.outFormatControl 0" );
	MGlobal::executeCommand( "setAttr defaultRenderGlobals.animation 1" );
	MGlobal::executeCommand( "setAttr defaultRenderGlobals.putFrameBeforeExt 1" ); // image name format
	MGlobal::executeCommand( "setAttr defaultRenderGlobals.periodInExt 1" );
	MGlobal::executeCommand( "setAttr defaultRenderGlobals.extensionPadding 1" ); // frame padding
	MGlobal::executeCommand( "setAttr defaultRenderGlobals.imageFormat 8" ); // jpg 8; bmp 20
	MGlobal::executeCommand( "setAttr defaultResolution.width 1280" ); // set image width
	MGlobal::executeCommand( "setAttr defaultResolution.height 720" ); // set image height

	std::string savePath = "";
	int folderNum = 1;

	do {
		savePath = Convenience::convertMStringToStdString(mllPath) + "/Render_Image_";
		Convenience::appendNumToStdString(savePath, folderNum);
		folderNum++;
	}
	while( Convenience::dirExists( savePath ) );

	_mkdir( savePath.c_str() );

	for ( unsigned int i = 1 ; i <= end_frame - start_frame + 1 ; ++i ) {
		std::stringstream ss;
		std::string frameCount = "";
		std::string frameIndex;

		if ( i < 10 ) {
			frameCount = "00";
		}
		else if ( i >= 10 && i < 100 ) {
			frameCount = "0";
		}
		else {
			frameCount = "";
		}

		Convenience::appendNumToStdString( frameCount, i );

		ss << start_frame + i - 1;
		ss >> frameIndex;
		ss.clear();

		std::string saveRender = ( std::string )"renderWindowSaveImageCallback " + '\"' + ( std::string )"renderView" + '\"' + " " + '\"' + savePath + "/Image_" + frameCount + '\"' + " "+ '\"'+ ( std::string )"JPEG"+ '\"';
		
		MString strSaveRender = ( char* )saveRender.c_str();
		MString strFrameIndex = ( char* )frameIndex.c_str();

		MGlobal::executeCommand( "currentTime " + strFrameIndex );

		MGlobal::executeCommand( "renderIntoNewWindow render" );
		MGlobal::executeCommand( strSaveRender );
	}

	return MStatus::kSuccess;
}