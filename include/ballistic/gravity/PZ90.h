#include "EarthModel.h"

namespace ball
{
	class PZ90 : public IEarth
	{
	public:
		size_t count() const override { return 36; }
		const std::vector<std::pair<double, double>>& harmonics() const { return _harmonics; }

		// a gravitational constant
		constexpr static inline double Mu() { return 398600.4418e9; }
		// an equatorial radius
		constexpr static inline double R() { return 6378136; }
		// a square of the eccentricity
		constexpr static inline double Ec2() { return 0.0066943662; }
		// a velocity of the rotation in rad/s
		constexpr static inline double W() { return 7.292115e-5; }
		// a flattening (a compression of Earth's ellipsoid)
		constexpr static inline double Fl() { return 1.0 / 298.2564151; }

	private:
		const static std::vector<std::pair<double, double>> _harmonics;
	};

}