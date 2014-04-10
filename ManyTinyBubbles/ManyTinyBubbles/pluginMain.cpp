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
#include "Convenience.h"


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
	std::string mll_filepath_std = Convenience::convertMStringToStdString( mll_filepath );
	std::string mel_filepath_std = mll_filepath_std.substr( 0, mll_filepath_std.find_last_of( "\\/" ) );
	MString mel_filepath = Convenience::convertStdStringToMString( mel_filepath_std );

	// load top-level menu script
	MString command = "source \"" + mel_filepath + "/gui.mel\";";
	MGlobal::executeCommand( command );

	// execute another script
	command = "source \"" + mel_filepath + "/create3dContainer.mel\";";
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

	// "many_tiny_bubbles_menu" is menu name defined in gui.mel
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