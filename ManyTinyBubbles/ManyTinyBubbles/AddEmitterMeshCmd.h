//#pragma once

#ifndef _AddEmitterMeshCmd
#define _AddEmitterMeshCmd

#include <maya/MPxCommand.h>


class AddEmitterMeshCmd : public MPxCommand
{

public:
	AddEmitterMeshCmd( void );
	virtual ~AddEmitterMeshCmd( void );

    static void* creator() { return new AddEmitterMeshCmd(); }

    MStatus doIt( const MArgList& args );
};

#endif