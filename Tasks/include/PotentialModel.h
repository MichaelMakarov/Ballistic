#pragma once
#include <utility>
#include <vector>

namespace ball
{
	namespace tasks
	{
		class IPotential
		{
		public:
			virtual unsigned int Count() const = 0;
			virtual const std::vector<std::pair<double, double>>& Harmonics() const = 0;
		};
	}
}