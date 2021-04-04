#pragma once
#include "general/Times.h"
#include <list>
#include <string>
#include <map>
#include <array>

namespace ball
{
	struct SpaceWeatherData
	{
		std::array<double, 8> Kp;
		std::array<double, 8> Ap;
		double Kpsum, Apavg, F10_7, F81;
	};

	std::map<general::time::Date, SpaceWeatherData> load_spaceweather_data(const std::string& filepath);

}