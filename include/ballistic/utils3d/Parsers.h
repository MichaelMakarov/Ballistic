#pragma once
#include <vector>
#include <string>
#include <fstream>
#include "Common.h"

namespace ball
{
	namespace graphics
	{
		std::ifstream& load_blender_model(std::ifstream& filein, std::vector<Vertex>& vertices);
	}
}