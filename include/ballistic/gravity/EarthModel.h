#pragma once
#include <utility>
#include <vector>

namespace ball
{
	/// <summary>
	/// Class represents Earth's gravity model and stores coefficients of geopotential model
	/// </summary>
	class IEarth
	{
	public:
		// a number of the potenital harmonics
		virtual size_t count() const = 0;
		// a list of the coefficicents (the potential harmonics)
		virtual const std::vector<std::pair<double, double>>& harmonics() const = 0;
	};
}