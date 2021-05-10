#pragma once
#include "general/Geometry.h"
#include "general/Times.h"
#include <utility>

namespace ball
{
	/// <summary>
	/// an acceleration caused in the point by single massive point
	/// </summary>
	/// <param name="point"> - a point to calculate an acceleration at</param>
	/// <param name="pos"> - a position of single massive point</param>
	/// <param name="mu"> - a gravitational parameter (multiplication of mass and gravitational const)</param>
	/// <returns>vector of dimension 3</returns>
	general::math::Vec3 acceleration_by_masspoint(
		const general::math::Vec3& point, 
		const general::math::Vec3& pos,
		const double mu);

	/// <summary>
	/// solar functionality representation
	/// </summary>
	class Sun
	{
	public:
		constexpr const static double Mu() { return 1.327124400189e20; }
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
		constexpr const static double Mu() { return 4.90486959e12; }
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