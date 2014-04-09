//
// Copyright (C) CIS660
// 
// File: pluginMain.cpp
//
// Author: Maya Plug-in Wizard 2.0
//

#include "ManyTinyBubblesNode.h"

#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>


////////////////////////////////////////////////////
// function prototypes
////////////////////////////////////////////////////

std::string convertMStringToStdString( MString mstring );
MString convertStdStringToMString( std::string std_string );


//
//	Description:
//		this method is called when the plug-in is loaded into Maya.  It 
//		registers all of the services that this plug-in provides with 
//		Maya.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
MStatus initializePlugin( MObject obj )
{
	MStatus   status;
	MFnPlugin plugin( obj, "CIS660", "2014", "Any" );


	////////////////////////////////////////////////////
	// register commands
	////////////////////////////////////////////////////

    //status = plugin.registerCommand( "CreateBubbleCmd", CreateBubbleCmd::creator );
    //if (!status) {
    //    status.perror("registerCommand");
    //    return status;
    //}

    //status = plugin.registerCommand( "RenderSequentialImagesCmd", RenderSequentialImagesCmd::creator );
    //if (!status) {
    //    status.perror("registerCommand");
    //    return status;
    //}


	////////////////////////////////////////////////////
	// register ManyTinyBubbles node
	////////////////////////////////////////////////////

	status = plugin.registerNode( "ManyTinyBubbles", ManyTinyBubbles::id, ManyTinyBubbles::creator,
								  ManyTinyBubbles::initialize );
	if ( !status ) {
		status.perror( "registerNode" );
		return status;
	}


	////////////////////////////////////////////////////
	// load menu options
	////////////////////////////////////////////////////

	// get .mll filepath
	MString mll_filepath = plugin.loadPath();

	// go up one level in path
	std::string mll_filepath_std = convertMStringToStdString( mll_filepath );
	std::string mel_filepath_std = mll_filepath_std.substr( 0, mll_filepath_std.find_last_of( "\\/" ) );
	MString mel_filepath = convertStdStringToMString( mel_filepath_std );

	// execute script
	MString command = "source \"" + mel_filepath + "/gui.mel\";";
	MGlobal::executeCommand( command );

	return status;
}


//
//	Description:
//		this method is called when the plug-in is unloaded from Maya. It 
//		deregisters all of the services that it was providing.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
MStatus uninitializePlugin( MObject obj )
{
	MStatus   status;
	MFnPlugin plugin( obj );


	////////////////////////////////////////////////////
	// remove top-level menu
	////////////////////////////////////////////////////

	// TODO: circumvent need for hardcoded menu name

	MGlobal::executeCommand( "if ( `menu -exists many_tiny_bubbles_menu` ) { deleteUI many_tiny_bubbles_menu; }" );


	////////////////////////////////////////////////////
	// deregister commands
	////////////////////////////////////////////////////

	//status = plugin.deregisterCommand( "CreateBubbleCmd" );
	//if ( !status ) {
	//	status.perror( "deregisterCommand" );
	//	return status;
	//}
	//
	//status = plugin.deregisterCommand( "RenderSequentialImagesCmd" );
	//if ( !status ) {
	//	status.perror( "deregisterCommand" );
	//	return status;
	//}


	////////////////////////////////////////////////////
	// deregister ManyTinyBubbles node
	////////////////////////////////////////////////////

	status = plugin.deregisterNode( ManyTinyBubbles::id );
	if ( !status ) {
		status.perror( "deregisterNode" );
		return status;
	}

	return status;
}



/*********** HELPER FUNCTIONS ***********/


////////////////////////////////////////////////////
// convertMStringToStdString
////////////////////////////////////////////////////
std::string convertMStringToStdString( MString mstring )
{
	std::string std_string = mstring.asChar();
	return std_string;
}


////////////////////////////////////////////////////
// convertStdStringToMString
////////////////////////////////////////////////////
MString convertStdStringToMString( std::string std_string )
{
	MString mstring = std_string.c_str();
	return mstring;
}