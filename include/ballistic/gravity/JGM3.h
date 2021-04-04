#pragma once
#include "EarthModel.h"

namespace ball
{
	class JGM3 : public IEarth
	{
	public:
		size_t count() const { return 50; }
		double Mu() const { return 0.3986004415E+15; }
		double R() const { return 0.6378136300E+07; }
		double Ec2() const { return 6.69437999014e-3; }
		double W() const { return 72.92115e-6; }
		double Fl() const { return 1.0 / 298.257223563; }
		const std::vector<std::pair<double, double>>& harmonics() const { return _harmonics; }

	private:
		const static std::vector<std::pair<double, double>> _harmonics;
	};
}