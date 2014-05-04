#pragma once

#ifndef RenderCmd_H_
#define RenderCmd_H_

#include <maya/MPxCommand.h>


class RenderCmd : public MPxCommand
{
public:
    RenderCmd( void );
    virtual ~RenderCmd( void );

    static void* creator() { return new RenderCmd(); }

    MStatus doIt( const MArgList& args );

	static MSyntax syntax();
};

#endif