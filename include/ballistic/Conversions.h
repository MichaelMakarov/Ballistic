#pragma once
#include "Structures.h"

namespace ball
{
	/// <summary>
	/// the height calculation from GCS position and the Earth's parameters
	/// </summary>
	/// <param name="pos">a point in GCS</param>
	/// <param name="rad">the Earth's equator radius</param>
	/// <param name="fl">the Earth's flatenning</param>
	/// <returns>a height above the Earth's ellipsoid</returns>
	double height_from_gcsposition(
		const general::math::Vec3& pos,
		const double rad, const double fl);
	// True anomaly calculation using the eccentric anomaly and the eccentricity.
	double trueanomaly_from_eccentric(const double E, const double e);
	// True anomaly calculation using the mean anomaly and the eccentricity.
	double trueanomaly_from_meananomaly(const double M, const double e);
	// The radius calculation using the semimajor axis, the true anomaly and the eccentricity.
	double radius_from_trueanomaly(const double a, const double v, const double e);
	// The rotation period calculating using a semimajor axis and gravity parameter
	double period_from_semimajoraxis(const double a, const double mu);
	// The semimajor axis calculating using a rotation period and gravity parameter
	double semimajoraxis_from_period(const double T, const double mu);

	/// <summary>
	/// Conversion from orthogonal coordinate system to spherical
	/// </summary>
	/// <param name="vec"> - a vector (x,y,z)</param>
	/// <returns>a vector (radius, latitude or inclination, longitude or ascension)</returns>
	general::math::Vec3 ort_to_sph(const general::math::Vec3& vec);
	/// <summary>
	/// Conversion the vector from spherical coordinate system to orthogonal
	/// </summary>
	/// <param name="vec"> - a vector  (radius, latitude or inclination, longitude or ascension)</param>
	/// <returns>a vector (x,y,z)</returns>
	general::math::Vec3 sph_to_ort(const general::math::Vec3& vec);
	/// <summary>
	/// Conversion from GCS to ACS
	/// </summary>
	/// <param name="vec"> - a vector (x,y,z)</param>
	/// <param name="t"> - sidereal time</param>
	/// <returns></returns>
	general::math::Vec3 ACS_to_GCS(const general::math::Vec3& vec, const double sidereal_time);
	/// <summary>
	/// Conversion from GCS to ACS
	/// </summary>
	/// <param name="vec"> - a vector (x,y,z)</param>
	/// <param name="t"> - sidereal time</param>
	/// <returns></returns>
	general::math::Vec3 GCS_to_ACS(const general::math::Vec3& vec, const double sidereal_time);
	/// <summary>
	/// conversion from ecliptic coordinate system to absolute
	/// </summary>
	/// <param name="vec"> - a vector (x, y, z)</param>
	/// <param name="e"> - ecliptic inclination</param>
	/// <returns>a vector (x, y, z) in ACS</returns>
	general::math::Vec3 ECS_to_ACS(const general::math::Vec3& vec, const double ecl_incl);
	/// <summary>
	/// conversion from absolute coordinate system to ecliptic
	/// </summary>
	/// <param name="vec"> - a vector (x, y, z)</param>
	/// <param name="e"> - ecliptic inclination</param>
	/// <returns>a vector (x, y, z) in ECS</returns>
	general::math::Vec3 ACS_to_ECS(const general::math::Vec3& vec, const double ecl_incl);
	/// <summary>
	/// calculating the julian centures since 2000
	/// </summary>
	/// <param name="jd"> - julian date related to mig=dnight</param>
	/// <returns></returns>
	inline constexpr double jd_to_jc2000(const general::time::JD& jd) noexcept {
		return (jd - general::time::JD2000 * general::time::SEC_PER_DAY).to_double() / 36525;
	}


	/// <summary>
	/// calculating orbital parameters from position and velocity
	/// </summary>
	/// <param name="pos">is a position vector in absolute cartesian coordinate system</param>
	/// <param name="vel">is a velocity vector in absolute cartesian coordinate system</param>
	/// <param name="mu">is a gravitational parameter</param>
	/// <returns>a struct of orbital parameters</returns>
	Oscul oscul_from_ACS(
		const general::math::Vec3& pos,
		const general::math::Vec3& vel,
		const double mu);
}