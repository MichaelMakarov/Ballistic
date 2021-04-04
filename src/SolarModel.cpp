#include "SolarModel.h"
#include "Conversions.h"
#include "general/Mathematics.h"

namespace ball
{
	general::math::Vec3 Sun::position_sphACS(const general::time::JD& jd)
	{
		using namespace general;
		const double t = jd_to_jc2000(jd);	// julian centures since 2000
		math::Vec3 coords;
		// a constant of the aberration
		constexpr double hi{ math::sec_to_rad(20.49552) };
		constexpr double ac{ 149597870691 };
		// solar average longitude
		const double L = math::sec_to_rad(1009677.85 + (100 * math::SEC_PER_ROUND + 2771.27 + 1.089 * t) * t);
		// solar perigee average longitude
		const double lc = math::sec_to_rad(1018578.046 + (6190.046 + (1.666 + 0.012 * t) * t) * t);
		// the Earth's orbit eccentricity
		const double e = 0.0167086342 - (0.000004203654 + (0.00000012673 + 0.00000000014 * t) * t) * t;
		// average ecliptic inclination
		const double eps = math::sec_to_rad(84381.448 - (46.815 + (0.00059 - 0.001813 * t) * t) * t);
		// ecliptic average longitude of lunar ascending node
		const double omega = math::sec_to_rad(450160.280 - (5 * math::SEC_PER_ROUND + 482890.539 - (7.455 + 0.008 * t) * t) * t);
		// long periodic nutation of the Earth
		const double psi = math::sec_to_rad(-17.1996 * std::sin(omega));
		// solar longitude
		const double longitude = L + 2 * e * std::sin(L - lc) + 1.25 * e * e * std::sin(2 * (L - lc));
		const double sinl{ std::sin(longitude) };
		const double cosl{ std::cos(longitude) };
		const double sine{ std::sin(eps) };
		const double cose{ std::cos(eps) };

		coords.Y = std::atan(sinl * sine / std::sqrt(cosl * cosl + sinl * sinl * cose * cose));
		coords.Z = std::atan(sinl / cosl * cose);

		if (coords.Z < 0.0) {
			if (coords.Y < 0.0) coords.Z += math::PI2;
			else coords.Z += math::PI;
		} else {
			if (coords.Y < 0.0) coords.Z += math::PI;
		}
		coords.Z += 0.061165 * psi - hi;
		coords.Y += hi * sine * cosl;
		coords.X = ac * (1 - e * (std::cos(L - lc) - e * 0.5 * (1 - std::cos(2 * (L - lc)))));

		return coords;
	}
	general::math::Vec3 Sun::position_ortACS(const general::time::JD& time)
	{
		return CS_sph_to_ort(position_sphACS(time));
	}
}