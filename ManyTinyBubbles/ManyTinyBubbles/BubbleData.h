#pragma once

#include "vec.h"
#include <vector>


class BubbleData
{
public:

	BubbleData( void );
	~BubbleData( void );

private:

	std::vector<std::vector<vec3>> bubblePosList;
	std::vector<std::vector<vec3>> bubbleVelList;
	std::vector<double> bubbleRadiusList;

};