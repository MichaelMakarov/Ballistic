#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
	/// <summary>
	/// class implements solar functionality
	/// </summary>
	class Sun
	{
	public:
		/// <summary>
		/// calculating the solar position in spherical ACS
		/// </summary>
		/// <param name="time">julian time reffered to midnight</param>
		/// <returns>vector (radius, inclination, ascension)</returns>
		static general::math::Vec3 position_sphACS(const general::time::JD& time);
		/// <summary>
		/// calculating the solar position in orthogonal ACS
		/// </summary>
		/// <param name="time">julian time reffered to midnight</param>
		/// <returns>vector (x, y, z)</returns>
		static general::math::Vec3 position_ortACS(const general::time::JD& time);
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