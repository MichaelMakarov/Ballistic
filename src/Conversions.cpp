#include "Conversions.h"
#include "general/Mathematics.h"

namespace ball
{
	double height_from_gcsposition(
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
	general::math::Vec3 ACS_to_GCS(const general::math::Vec3& vec, const double sidereal_time)
	{
		const double sint = std::sin(sidereal_time), cost = std::cos(sidereal_time);
		return general::math::Vec3{
			vec[0] * cost + vec[1] * sint,
			vec[1] * cost - vec[0] * sint,
			vec[2] };
	}
	general::math::Vec3 GCS_to_ACS(const general::math::Vec3& vec, const double sidereal_time)
	{
		const double sint = std::sin(sidereal_time), cost = std::cos(sidereal_time);
		return general::math::Vec3{
			vec[0] * cost - vec[1] * sint,
			vec[1] * cost + vec[0] * sint,
			vec[2] };
	}
	general::math::Vec3 ECS_to_ACS(const general::math::Vec3& vec, const double e)
	{
		const double sine = std::sin(e), cose = std::cos(e);
		return general::math::Vec3{
			vec[0],
			vec[1] * cose - vec[2] * sine,
			vec[1] * sine + vec[2] * cose
		};
	}
	general::math::Vec3 ACS_to_ECS(const general::math::Vec3& vec, const double e)
	{
		const double sine = std::sin(e), cose = std::cos(e);
		return general::math::Vec3{
			vec[0],
			vec[2] * sine + vec[1] * cose,
			vec[2] * cose - vec[1] * sine
		};
	}

	Oscul oscul_from_ACS(
		const general::math::Vec3& pos,
		const general::math::Vec3& vel,
		const double mu) 
	{
		using namespace general::math;
		Oscul osc;
		const Vec3 h = cross(pos, vel);
		const auto w = normalize(h);
		const double p{ h * h / mu };
		const double r{ pos.length() };
		osc.inclination = std::atan(std::sqrt(w[0] * w[0] + w[1] * w[1]) / w[2]);
		osc.ascendnode = std::atan2(w[0], -w[1]);
		osc.semiaxis = 1 / (2 / pos.length() - vel * vel / mu);
		osc.eccentricity = std::sqrt(1 - p / osc.semiaxis);
		osc.ecanomaly = std::atan((pos * vel) / std::sqrt(osc.semiaxis * mu) / (1 - r / osc.semiaxis));
		osc.meananomaly = osc.ecanomaly - osc.eccentricity * std::sin(osc.ecanomaly);
		osc.trueanomaly = std::atan2(
			std::sqrt(1 - osc.eccentricity * osc.eccentricity) * std::sin(osc.ecanomaly),
			std::cos(osc.ecanomaly) - osc.eccentricity);
		osc.latitudearg = std::atan2(pos[2], pos[1] * w[0] - pos[0] * w[1]);
		osc.periapsis = osc.trueanomaly - osc.latitudearg;
		return osc;
	}

}