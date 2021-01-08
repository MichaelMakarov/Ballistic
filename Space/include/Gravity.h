#pragma once
#include <utility>
#include <vector>

namespace ball
{
	namespace space
	{
		class IGravity
		{
		public:
			// a gravitational constant
			virtual double Mu() const = 0;
			// an equatorial radius
			virtual double R() const = 0;
			// a square of the eccentricity
			virtual double Ec2() const = 0;
			// a velocity of the rotation in rad/s
			virtual double W() const = 0;
			// a flattening (a compression of Earth's ellipsoid)
			virtual double Fl() const = 0;
			// a number of the potenital harmonics
			virtual size_t Count() const = 0;
			// a list of the coefficicents (the potential harmonics)
			virtual const std::vector<std::pair<double, double>>& Harmonics() const = 0;
		};
	}
}