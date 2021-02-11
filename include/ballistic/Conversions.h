#pragma once
#include "general/GeneralConstants.h"
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
	namespace space
	{
		// Conversion julian date to julian centures since 2000/1/1 12:00:00.
		// jd - julian date related to midnight.
		inline double JD2000_date_to_centures(const general::time::JD& jd)
		{
			return ((jd - general::time::JD2000).to_double()) / 36525;
		}

		// The height calculation from GCS position and geo parameters.
		// rad - a radius of Earth's equator.
		// fl - the flatening of Earth.
		inline double GCS_height_from_position(
			const double x, const double y, const double z,
			const double rad, const double fl)
		{
			double dist = std::sqrt(x * x + y * y + z * z);
			return dist - rad * (1.0 - fl * z * z / dist / dist);
		}
		// The height calculation from GCS position and geo parameters.
		// rad - a radius of Earth's equator.
		// fl - the flatening of Earth.
		inline double GCS_height_from_position(
			const general::math::Vec3& pos,
			const double rad, const double fl)
		{
			double dist = pos.length();
			return dist - rad * (1.0 - fl * pos.Z * pos.Z / dist / dist);
		}

		// True anomaly calculation using the eccentric anomaly and the eccentricity.
		inline double trueanomaly_from_eccentric(const double E, const double e)
		{
			double cosv = (std::cos(E) - e) / (1 - e * std::cos(E)),
				sinv = std::sqrt(1 - e * e) * std::sin(E) / (1 - e * std::cos(E));
			return std::atan2(sinv, cosv);
		}
		// True anomaly calculation using the mean anomaly and the eccentricity.
		inline double trueanomaly_from_meananomaly(const double M, const double e)
		{
			return M + e * ((2 - 0.25 * e * e) * std::sin(M) +
				e * (1.25 * std::sin(2 * M) +
					13 / 12 * e * std::sin(3 * M)));
		}

		// The radius calculation using the semimajor axis, the true anomaly and the eccentricity.
		inline double radius_from_trueanomaly(const double a, const double v, const double e)
		{
			return a * (1 - e * e) / (1 + e * std::cos(v));
		}

		// The rotation period calculating using a semimajor axis and gravity parameter
		inline double period_from_semimajoraxis(const double a, const double mu) 
		{
			return 2 * general::math::PI * std::sqrt(a * a * a / mu);
		}
		// The semimajor axis calculating using a rotation period and gravity parameter
		inline double semimajoraxis_from_period(const double T, const double mu)
		{
			return std::pow(mu * T * T / (4 * general::math::PI * general::math::PI), 1.0 / 3);
		}

		// Convertion the vector from orthogonal coordinate system to spherical
		inline general::math::Vec3 CS_ortho_to_spher(const general::math::Vec3& xyz)
		{
			return general::math::Vec3(
				xyz.length(),
				std::atan2(xyz.Z, std::sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y)),
				std::atan2(xyz.Y, xyz.X));
		}
		// Convertion the vector from spherical coordinate system to orthogonal
		// rbl = (R, B, L) vector
		inline general::math::Vec3 CS_spher_to_ortho(const general::math::Vec3& rbl)
		{
			const double cosB = std::cos(rbl.Y);
			return general::math::Vec3(
				rbl.X * std::cos(rbl.Z) * cosB,
				rbl.X * std::sin(rbl.Z) * cosB,
				rbl.X * std::sin(rbl.Y));
		}
		// Convertion from absolute to Grinweech orthogonal coordinate system.
		// t - astro time
		inline general::math::Vec3 ACS_to_GCS(const general::math::Vec3& coords, const double t)
		{
			const double sint{ std::sin(t) }, cost{ std::cos(t) };
			return general::math::Vec3(
				coords.X * cost + coords.Y * sint,
				coords.Y * cost - coords.X * sint,
				coords.Z);
		}
		// Convertion from Grinweech to absolute orthogonal coordinate system.
		// t - astro time
		inline general::math::Vec3 GCS_to_ACS(const general::math::Vec3& coords, const double t)
		{
			const double sint{ std::sin(t) }, cost{ std::cos(t) };
			return general::math::Vec3(
				coords.X * cost - coords.Y * sint,
				coords.Y * cost + coords.X * sint,
				coords.Z);
		}
	}
}