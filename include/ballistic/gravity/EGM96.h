﻿#pragma once
#include "EarthModel.h"

namespace ball
{
	class EGM96 : public IEarth
	{
	public:
		size_t count() const override { return 50; }
		double Mu() const override { return 0.3986004415E+15; }
		double R() const override { return 0.6378136300E+07; }
		double Ec2() const override { return 6.69437999014e-3; }
		double W() const override { return 72.92115e-6; }
		double Fl() const override { return 1.0 / 298.257223563; }
		const std::vector<std::pair<double, double>>& harmonics() const { return _harmonics; }

	private:
		const static std::vector<std::pair<double, double>> _harmonics;
	};
}