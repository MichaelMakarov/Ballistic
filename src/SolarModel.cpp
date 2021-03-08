#include "SolarModel.h"
#include "general/Mathematics.h"
#include "Conversions.h"

namespace ball
{
	general::math::Vec3 SphACS_solar_position(const double T)
	{
		using namespace general;
		math::Vec3 coords;
		// a constant of the aberration
		const double hi{ math::sec_to_rad(20.49552) };
		const double ac{ 1.4959787e11 };
		// solar average longitude
		const double L = math::sec_to_rad(1009677.85 + (100 * math::SEC_PER_ROUND + 2771.27 + 1.089 * T) * T);
		// solar perigee average longitude
		const double lc = math::sec_to_rad(1018578.046 + (6190.046 + (1.666 + 0.012 * T) * T) * T);
		// Earth orbit eccentricity
		const double e = 0.0167086342 - (0.000004203654 + (0.00000012673 + 0.00000000014 * T) * T) * T;
		// ecliptic average inclination
		const double eps = math::sec_to_rad(84381.448 - (46.815 + (0.00059 - 0.001813 * T) * T) * T);
		// ecliptic average longitude of lunar ascending node
		const double omega = math::sec_to_rad(450160.280 - (5 * math::SEC_PER_ROUND + 482890.539 - (7.455 + 0.008 * T) * T) * T);
		// long periodic nutation of the Earth
		const double psi = math::sec_to_rad(-17.1996 * std::sin(omega));
		// solar longitude
		const double longitude = L + 2 * e * std::sin(L - lc) + 1.25 * e * e * std::sin(2 * (L - lc));
		const double sinL{ std::sin(longitude) };
		const double cosL{ std::cos(longitude) };
		const double sinE{ std::sin(eps) };
		const double cosE{ std::cos(eps) };

		coords.Y = std::atan(sinL * sinE / std::sqrt(cosL * cosL + sinL * sinL * cosE * cosE));
		coords.Z = std::atan(sinL / cosL * cosE);

		if (coords.Z < 0.0) {
			if (coords.Y < 0.0) coords.Z += math::PI2;
			else coords.Z += math::PI;
		}
		else {
			if (coords.Y < 0.0) coords.Z += math::PI;
		}
		coords.Z += psi - hi;
		coords.Y += hi * sinE * cosL;
		coords.X = ac * (1 - eps * (std::cos(L - lc) + eps * 0.5 * (1 - std::cos(2 * (L - lc)))));

		return coords;
	}
	general::math::Vec3 ACS_solar_position(const double T)
	{
		return CS_spher_to_ortho(SphACS_solar_position(T));
	}
}