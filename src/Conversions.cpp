#include "Conversions.h"
#include "general/Mathematics.h"

namespace ball
{
	double GCS_height_from_position(
		const general::math::Vec3& pos,
		const double rad, const double fl)
	{
		const double dist = pos.length();
		return dist - rad * (1.0 - fl * pos[2] * pos[2] / dist / dist);
	}
	double trueanomaly_from_eccentric(const double E, const double e)
	{
		const double cosv = (std::cos(E) - e) / (1 - e * std::cos(E)),
			sinv = std::sqrt(1 - e * e) * std::sin(E) / (1 - e * std::cos(E));
		return std::atan2(sinv, cosv);
	}
	double trueanomaly_from_meananomaly(const double M, const double e)
	{
		return M + e * ((2 - 0.25 * e * e) * std::sin(M) +
			e * (1.25 * std::sin(2 * M) +
				13 / 12 * e * std::sin(3 * M)));
	}
	double radius_from_trueanomaly(const double a, const double v, const double e)
	{
		return a * (1 - e * e) / (1 + e * std::cos(v));
	}
	double period_from_semimajoraxis(const double a, const double mu)
	{
		return 2 * general::math::PI * std::sqrt(a * a * a / mu);
	}
	double semimajoraxis_from_period(const double T, const double mu)
	{
		return std::pow(mu * T * T / (4 * general::math::PI * general::math::PI), 1.0 / 3);
	}

	general::math::Vec3 ort_to_sph(const general::math::Vec3& vec)
	{
		return general::math::Vec3({
			vec.length(),
			std::atan2(vec[2], std::sqrt(vec[0] * vec[0] + vec[1] * vec[1])),
			std::atan2(vec[1], vec[0]) });
	}
	general::math::Vec3 sph_to_ort(const general::math::Vec3& vec)
	{
		const double cosB = std::cos(vec[1]);
		return general::math::Vec3({
			vec[0] * std::cos(vec[2]) * cosB,
			vec[0] * std::sin(vec[2]) * cosB,
			vec[0] * std::sin(vec[1]) });
	}
	general::math::Vec3 ACS_to_GCS(const general::math::Vec3& vec, const double t)
	{
		const double sint{ std::sin(t) }, cost{ std::cos(t) };
		return general::math::Vec3({
			vec[0] * cost + vec[1] * sint,
			vec[1] * cost - vec[0] * sint,
			vec[2] });
	}
	general::math::Vec3 GCS_to_ACS(const general::math::Vec3& vec, const double t)
	{
		const double sint{ std::sin(t) }, cost{ std::cos(t) };
		return general::math::Vec3({
			vec[0] * cost - vec[1] * sint,
			vec[1] * cost + vec[0] * sint,
			vec[2] });
	}
	general::math::Vec3 ECS_to_ACS(const general::math::Vec3& vec, const double e)
	{
		const double sine{ std::sin(e) }, cose{ std::cos(e) };
		return general::math::Vec3({
			vec[0],
			vec[1] * cose - vec[2] * sine,
			vec[1] * sine + vec[2] * cose }
		);
	}
	general::math::Vec3 ACS_to_ECS(const general::math::Vec3& vec, const double e)
	{
		const double sine{ std::sin(e) }, cose{ std::cos(e) };
		return general::math::Vec3({
			vec[0],
			vec[2] * sine + vec[1] * cose,
			vec[2] * cose - vec[1] * sine
		});
	}
	double jd_to_jc2000(const general::time::JD& jd)
	{
		return (jd - (general::time::JD2000 + 0.5)).to_double() / 36525;
	}

	double sidereal_time(const double m, const double t)
	{
		using namespace general::math;
		return rad_to_2pi(1.7533685592 + 6.2831853072 * m +	t * (0.0172027918051 * 36525 + t * (6.7707139e-6 - 4.50876e-10 * t)));
	}
	double sidereal_time_avr(const general::time::JD& jd, const double timezone)
	{
		return sidereal_time(jd.T() - timezone / 24.0, jd_to_jc2000(jd - timezone / 24.0));
	}
	double sidereal_time_true(const general::time::JD& jd, const double timezone)
	{
		using namespace general::math;
		double t{ jd_to_jc2000(jd - timezone / 24.0) };
		// ecliptic average longitude of lunar ascending node
		const double omega = sec_to_rad(450160.280 - (5 * SEC_PER_ROUND + 482890.539 - (7.455 + 0.008 * t) * t) * t);
		// the Earth's nutation in ascension
		const double nut = 0.061165 * sec_to_rad(-17.1996 * std::sin(omega));
		return rad_to_2pi(sidereal_time_avr(jd, timezone) + nut);
	}

}