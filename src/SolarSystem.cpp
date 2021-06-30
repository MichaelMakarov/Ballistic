#include "SolarSystem.h"

namespace ball
{
	std::pair<general::math::Vec3, general::math::Vec3> Sun::positionACS(const general::time::JD& jd)
	{
		using namespace general::math;
		const double t = jd_to_jc2000(jd);	// julian centures since 2000
		Vec3 coords;
		// a constant of the aberration
		constexpr double hi{ sec_to_rad(20.49552) };
		constexpr double ac{ 149597870691 };
		// solar average longitude
		const double L = sec_to_rad(1009677.85 + (100 * SEC_PER_ROUND + 2771.27 + 1.089 * t) * t);
		// solar perigee average longitude
		const double lc = sec_to_rad(1018578.046 + (6190.046 + (1.666 + 0.012 * t) * t) * t);
		// the Earth's orbit eccentricity
		const double e = 0.0167086342 - (0.000004203654 + (0.00000012673 + 0.00000000014 * t) * t) * t;
		// average ecliptic inclination
		const double eps = sec_to_rad(84381.448 - (46.815 + (0.00059 - 0.001813 * t) * t) * t);
		// ecliptic average longitude of lunar ascending node
		const double omega = sec_to_rad(450160.280 - (5 * SEC_PER_ROUND + 482890.539 - (7.455 + 0.008 * t) * t) * t);
		// long periodic nutation of the Earth
		const double psi = sec_to_rad(-17.1996 * std::sin(omega));
		// solar longitude
		const double longitude = L + 2 * e * std::sin(L - lc) + 1.25 * e * e * std::sin(2 * (L - lc));
		const double sinl{ std::sin(longitude) };
		const double cosl{ std::cos(longitude) };
		const double sine{ std::sin(eps) };
		const double cose{ std::cos(eps) };

		coords[1] = std::atan(sinl * sine / std::sqrt(cosl * cosl + sinl * sinl * cose * cose));
		coords[2] = std::atan(sinl / cosl * cose);

		if (coords[2] < 0.0) {
			if (coords[1] < 0.0) coords[2] += PI2;
			else coords[2] += PI;
		} else {
			if (coords[1] < 0.0) coords[2] += PI;
		}
		coords[2] += 0.061165 * psi - hi;
		coords[1] += hi * sine * cosl;
		coords[0] = ac * (1 - e * (std::cos(L - lc) - e * 0.5 * (1 - std::cos(2 * (L - lc)))));

		return std::make_pair(coords, sph_to_ort(coords));
	}
	general::math::Vec3 Sun::acceleration(const general::math::Vec3& point, const general::time::JD& time)
	{
		auto [tmp, pos] = positionACS(time);
		pos = ACS_to_GCS(pos, sidereal_time_true(time));
		return acceleration_by_masspoint(point, pos, Mu());
	}

	std::pair<general::math::Vec3, general::math::Vec3> Moon::positionACS(const general::time::JD& jd)
	{
		using namespace general::math;
		const double t = jd_to_jc2000(jd);	// julian centures since 2000
		// radius of Earth's equator
		const double r{ 6378136 };
		// average lunar anomaly
		const double la = sec_to_rad((485866.733 + (1325 * SEC_PER_ROUND + 715922.633 + (31.31 + 0.064 * t) * t) * t));
		// solar average anomaly
		const double sa = sec_to_rad((1287099.804 + (99 * SEC_PER_ROUND + 1292581.224 - (0.577 + 0.012 * t) * t) * t));
		// average arg of lunar latitude
		const double f = sec_to_rad((335778.877 + (1342 * SEC_PER_ROUND + 295263.137 - (13.257 - 0.011 * t) * t) * t));
		// average alongation (a difference between solar and lunar longitudes)
		const double d = sec_to_rad((1072261.307 + (1236 * SEC_PER_ROUND + 1105601.328 - (6.891 - 0.019 * t) * t) * t));
		// lunar ecliptic latitude
		double latitude = sec_to_rad(18461.48 * std::sin(f) +
			1010.18 * std::sin(la + f) - 999.69 * std::sin(f - la) -
			623.65 * std::sin(f - 2 * d) + 199.48 * std::sin(f + 2 * d - la) -
			166.57 * std::sin(la + f - 2 * d) +
			117.26 * std::sin(f + 2 * d) + 61.91 * std::sin(2 * la + f) -
			33.35 * std::sin(f - 2 * d - la) - 31.76 * std::sin(f - 2 * la) -
			29.68 * std::sin(sa + f - 2 * d) + 15.125 * std::sin(la + f + 2 * d) -
			15.56 * std::sin(2 * (la - d) + f));
		// lunar ecliptic longitude
		double longitude = sec_to_rad(785939.157 +
			(1336 * SEC_PER_ROUND + 1108372.598 + (5.802 + 0.019 * t) * t) * t +
			22639.5 * std::sin(la) - 4586.42 * std::sin(la - 2 * d) + 2369.9 * std::sin(2 * d) +
			769.01 * std::sin(2 * la) - 668.11 * std::sin(sa) - 411.6 * std::sin(2 * f) -
			211.65 * std::sin(2 * (la - d)) -
			205.96 * std::sin(la + sa - 2 * d) + 191.95 * std::sin(la + 2 * d) -
			165.14 * std::sin(sa - 2 * d) +
			147.69 * std::sin(la - sa) - 125.15 * std::sin(d) - 109.66 * std::sin(la + sa) -
			55.17 * std::sin(2 * (f - d)) - 45.1 * std::sin(sa + 2 * f) + 39.53 * std::sin(la - 2 * f) -
			38.42 * std::sin(la - 4 * d) + 36.12 * std::sin(3 * la) - 30.77 * std::sin(2 * la - 4 * d) +
			28.47 * std::sin(la - sa - 2 * d) - 24.42 * std::sin(sa + 2 * d) + 18.6 * std::sin(la - d) +
			18.02 * std::sin(sa - d));
		longitude = rad_to_2pi(longitude);
		// lunar paralax
		double paralax = sec_to_rad(3422.7 + 186.539 * std::cos(la) + 34.311 * std::cos(la - 2 * d) +
			28.233 * std::cos(2 * d) + 10.165 * std::cos(2 * la) + 3.086 * std::cos(la + 2 * d) +
			1.92 * std::cos(sa - 2 * d) + 1.445 * std::cos(la + sa - 2 * d) + 1.154 * std::cos(la - sa) -
			0.975 * std::cos(d) - 0.95 * std::cos(la + sa) - 0.713 * std::cos(la - 2 * f) +
			0.6215 * std::cos(3 * la) + 0.601 * std::cos(la - 4 * d));
		double radius = r / paralax;
		Vec3 coords = ECS_to_ACS(
			sph_to_ort(Vec3({ radius, latitude, longitude })), 
			sec_to_rad(84381.448 - (46.815 + (0.00059 - 0.001813 * t) * t) * t));
		return std::make_pair(ort_to_sph(coords), coords);
	}
	general::math::Vec3 Moon::acceleration(const general::math::Vec3& point, const general::time::JD& time)
	{
		auto [tmp, pos] = positionACS(time);
		pos = ACS_to_GCS(pos, sidereal_time_true(time));
		return acceleration_by_masspoint(point, pos, Mu());
	}
	
	general::math::Vec3 acceleration_by_masspoint(const general::math::Vec3& point, const general::math::Vec3& pos, const double mu)
	{
		auto diff = pos - point;
		return mu * (diff / std::pow(diff.length(), 3) - pos / std::pow(pos.length(), 3));
	}

	general::math::Vec3 accelerationdiff_by_masspoint(const general::math::Vec3& point, const general::math::Vec3& pos, const double mu)
	{
		using namespace general::math;
		auto diff{ pos - point };
		const double dist = diff.length();
		diff[0] *= diff[0];
		diff[1] *= diff[1];
		diff[2] *= diff[2];
		return mu / std::pow(dist, 3) * (diff * 3 / point.length() / dist - Vec3::ones());
	}

	double sidereal_time_true(const general::time::JD& jd) noexcept
	{
		using namespace general::math;
		double jc = jd_to_jc2000(jd);
		// ecliptic inclination
		double e = Sun::ecliptic_mean_incl(jc);
		// lunar mean anomaly
		double la = Moon::ecl_mean_anomaly(jc);
		//solar mean anomaly
		double sa = Sun::ecl_mean_anomaly(jc);
		// lunar mean argument of latitude
		double f = Moon::ecl_mean_latarg(jc);
		// difference between lunar and solar longitudes
		double d = Sun::ecl_delta_lslong(jc);
		// ecliptic mean longitude of lunar ascending node
		double o = Moon::ecl_mean_ascnode_long(jc);
		// the Earth's nutation in ascension
		double nut = -0.83386e-4 * std::sin(o) + 0.9997e-6 * std::sin(2 * o) + 0.6913e-6 * std::sin(sa) -
			0.63932e-5 * std::sin(2 * (f - d + o)) - 0.11024e-5 * std::sin(2 * (f + o));
		return rad_to_2pi(sidereal_time_mean(jd) + nut * std::cos(e));
	}
	
}