#pragma once
#include "Constants.h"
#include "RBL.h"
#include "XYZ.h"
#include <cmath>

namespace ball
{
	namespace tasks
	{
		// True anomaly calculation using the eccentric anomaly and the eccentricity.
		inline double TrueAnomalyFromEccentric(const double E, const double e)
		{
			double cosv = (std::cos(E) - e) / (1 - e * std::cos(E)),
				sinv = std::sqrt(1 - e * e) * std::sin(E) / (1 - e * std::cos(E));
			return std::atan2(sinv, cosv);
		}
		// True anomaly calculation using the mean anomaly and the eccentricity.
		double TrueAnomalyFromMeanAnomaly(const double M, const double e)
		{
			return M + e * ((2 - 0.25 * e * e) * std::sin(M) +
				e * (1.25 * std::sin(2 * M) +
					13 / 12 * e * std::sin(3 * M)));
		}

		// The radius calculation using the semimajor axis, the true anomaly and the eccentricity.
		double RadiusFromTrueAnomaly(const double a, const double v, const double e)
		{
			return a * (1 - e * e) / (1 + e * std::cos(v));
		}

		// The rotation period calculating using a semimajor axis and gravity parameter
		inline double PeriodFromSemimajorAxis(const double a, const double mu) 
		{
			return 2 * math::PI * std::sqrt(a * a * a / mu);
		}
		// The semimajor axis calculating using a rotation period and gravity parameter
		inline double SemimajorAxisFromPeriod(const double T, const double mu)
		{
			return std::pow(mu * T * T / (4 * math::PI * math::PI), 1.0 / 3);
		}

		// Convertion the vector from orthogonal geocentral coordinate system to spherical
		inline geometry::RBL GCS_OrthoToSpher(const geometry::XYZ& xyz)
		{
			return geometry::RBL(
				xyz.Length(),
				std::atan2(xyz.Z, std::sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y)),
				std::atan2(xyz.Y, xyz.X));
		}
		// Convertion the vector from spherical geocentral coordinate system to orthogonal
		inline geometry::XYZ GCS_SpherToOrtho(const geometry::RBL& rbl)
		{
			double cosB = std::cos(rbl.B);
			return geometry::XYZ(
				rbl.R * std::cos(rbl.L) * cosB,
				rbl.R * std::sin(rbl.L) * cosB,
				rbl.R * std::sin(rbl.B));
		}
	}
}