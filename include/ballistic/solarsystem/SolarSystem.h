#pragma once
#include "Conversions.h"
#include <utility>

namespace ball
{
	/// <summary>
	/// an acceleration caused by single massive point
	/// </summary>
	/// <param name="point"> - a point to calculate an acceleration at</param>
	/// <param name="masspoint"> - a position of single massive point</param>
	/// <param name="mu"> - a gravitational parameter (multiplication of mass and gravitational const)</param>
	/// <returns>vector (dU/dx, dU/dy, dU/dz)</returns>
	general::math::Vec3 acceleration_by_masspoint(
		const general::math::Vec3& point, 
		const general::math::Vec3& masspoint,
		const double mu);
	/// <summary>
	/// a vector of second partial acceleration of the single massive point potential
	/// </summary>
	/// <param name="point">is a point to calculate at</param>
	/// <param name="masspoint">is a position of a single massive point</param>
	/// <param name="mu">is gravitational parameter</param>
	/// <returns>vector (ddU/dx/dx, ddU/dy/dy, ddU/dz/dz)<returns>
	general::math::Vec3 accelerationdiff_by_masspoint(
		const general::math::Vec3& point,
		const general::math::Vec3& masspoint,
		const double mu
	);

	/// <summary>
	/// average sidereal time according RD50
	/// </summary>
	/// <param name="jd">is a julian date refered to midnight</param>
	/// <returns>sidereal time in rad</returns>
	inline constexpr double sidereal_time_mean(const general::time::JD& jd) noexcept {
		using namespace general::math;
		double jc = jd_to_jc2000(jd);
		return rad_to_2pi(1.7533685592 + 6.2831853072 * jd.T() + jc * (0.0172027918051 * 36525 + jc * (6.7707139e-6 - 4.50876e-10 * jc)));
	}

	double sidereal_time_true(const general::time::JD& jd) noexcept;

	/// <summary>
	/// solar functionality representation
	/// </summary>
	class Sun
	{
	public:
		constexpr static double Mu() { return 1.327124400189e20; }
		/// <summary>
		/// eclipatic solar mean anomaly
		/// </summary>
		/// <param name="jc">- julian centures refered 2000 epoch</param>
		/// <returns></returns>
		constexpr static double ecl_mean_anomaly(const double jc) {
			return 6.24003594 + (628.30195602 - (2.7974e-6 + 5.82e-8 * jc) * jc) * jc;
		}
		/// <summary>
		/// ecliptic mean inclination
		/// </summary>
		/// <param name="jc">- julian centures refered 2000 epoch</param>
		/// <returns></returns>
		constexpr static double ecliptic_mean_incl(const double jc) {
			return 0.4090928042 - (0.2269655e-3 + (0.29e-8 - 0.88e-8 * jc) * jc) * jc;
		}
		/// <summary>
		/// difference between lunar and solar ecliptic longitudes
		/// </summary>
		/// <param name="jc">- julian centures refered 2000 epoch</param>
		/// <returns></returns>
		constexpr static double ecl_delta_lslong(const double jc) {
			return 5.19846951 + (7771.37714617 - (3.34085e-5 - 9.21e-8 * jc) * jc) * jc;
		}
		/// <summary>
		/// solar position in spherical and orthogonal ACS
		/// </summary>
		/// <param name="time"> - julian time reffered to midnight</param>
		/// <returns>pair of vector (radius, inclination, ascension) and vector (x, y, z)</returns>
		static std::pair<general::math::Vec3, general::math::Vec3> positionACS(const general::time::JD& time);
		/// <summary>
		/// an acceleration caused by Sun
		/// </summary>
		/// <param name="point"> - a point to calculate an acceleration at</param>
		/// <param name="time"> - julian time reffered to midnight</param>
		/// <returns>vector in GCS</returns>
		static general::math::Vec3 acceleration(const general::math::Vec3& point, const general::time::JD& time);
	};
	/// <summary>
	/// lunar functionality representation
	/// </summary>
	class Moon
	{
	public:
		/// <summary>
		/// gravitationsl parameter
		/// </summary>
		constexpr static double Mu() { return 4.90486959e12; }
		/// <summary>
		/// ecliptic luanr mean anomaly
		/// </summary>
		/// <param name="jc">- julian centures refered 2000 epoch</param>
		/// <returns></returns>
		constexpr static double ecl_mean_anomaly(const double jc) {
			return 2.355548393 + (8328.69142288 + (1.517952e-1 + 3.103e-7 * jc) * jc) * jc;
		}
		/// <summary>
		/// ecliptic lunar argument of latitude
		/// </summary>
		/// <param name="jc">- julian centures refered 2000 epoch</param>
		/// <returns></returns>
		constexpr static double ecl_mean_latarg(const double jc) {
			return 1.62790193 + (8433.46615831 - (6.42717e-5 - 5.33e-8 * jc) * jc) * jc;
		}
		/// <summary>
		/// ecliptic mean longitude of lunar ascending node
		/// </summary>
		/// <param name="jc">- julian centures refered 2000 epoch</param>
		/// <returns></returns>
		constexpr static double ecl_mean_ascnode_long(const double jc) {
			return 2.182438624 - (33.757045936 - (3.61429e-5 + 3.88e-8 * jc) * jc) * jc;
		}
		/// <summary>
		/// difference between lunar and solar ecliptic longitudes
		/// </summary>
		/// <param name="jc">- julian centures refered 2000 epoch</param>
		/// <returns></returns>
		constexpr static double ecl_delta_lslong(const double jc) {
			return 5.19846951 + (7771.37714617 - (3.34085e-5 - 9.21e-8 * jc) * jc) * jc;
		}
		/// <summary>
		/// lunar position in spherical and orthogonal ACS
		/// </summary>
		/// <param name="time"> - julian time reffered to midnight</param>
		/// <returns>pair of vector (radius, inclination, ascension) and vector (x, y, z)</returns>
		static std::pair<general::math::Vec3, general::math::Vec3> positionACS(const general::time::JD& time);
		/// <summary>
		/// an acceleration caused by Moon
		/// </summary>
		/// <param name="point"> - a point to calculate an acceleration at</param>
		/// <param name="time">julian time reffered to midnight</param>
		/// <returns>vector in GCS</returns>
		static general::math::Vec3 acceleration(const general::math::Vec3& point, const general::time::JD& time);
	};

	// Solar gravity model
	template<class SungravityModel>
	class ISungravity
	{
	public:
		general::math::Vec3 acceleration(
			const general::math::Vec3& solarpos,
			const general::math::Vec3& pointpos) const
		{
			return static_cast<const SungravityModel>(this)->acceleration(solarpos, pointpos);
		}
	};

	// Solar pressure model
	template<class SunlightModel, class ... Args>
	class ISunlight
	{
	public:
		general::math::Vec3 acceleration(
			const general::math::Vec3& solarpos,
			const general::math::Vec3& pointpos,
			const Args ... args) const
		{
			return static_cast<const SunlightModel>(this)->acceleration(solarpos, pointpos, args ...);
		}
	};
}