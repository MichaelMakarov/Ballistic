#include "EGM08.h"
#include <fstream>

namespace ball
{
	EGM08::EGM08(const char* filepath)
	{
		auto fin = std::ifstream(filepath);
		load_harmonics<2190>(fin, _harmonics, _count);
		fin.close();
	}

	size_t EGM08::_count = 50;

	std::vector<std::pair<double, double>> EGM08::_harmonics = std::vector<std::pair<double, double>>();
}