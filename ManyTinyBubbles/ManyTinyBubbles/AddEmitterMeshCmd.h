#pragma once

#ifndef AddEmitterMeshCmd_H_
#define AddEmitterMeshCmd_H_

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