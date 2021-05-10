#pragma once
#include <utility>
#include <vector>
#include <istream>

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
	/// <summary>
	/// loading harmonics from stream
	/// </summary>
	/// <param name="is">input stream</param>
	/// <param name="list">a list of harmonics</param>
	/// <param name="count">a number of harmonics</param>
	/// <param name="number">an expected number of harmonics</param>
	template<size_t number> void load_harmonics(
		std::istream& is, 
		std::vector<std::pair<double, double>>& list, 
		size_t& count)
	{
		if (!is) return;
		const size_t size = ((number + 2) * (number + 1)) / 2;
		list.resize(size);
		count = 0;
		while (!is.eof() && count < size) {
			is >> list[count].first >> list[count].second;
			++count;
		}
		count = number;
	}
}