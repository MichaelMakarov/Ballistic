#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
	// Calculating the solar position in spherical absolute coordinate system
	general::math::Vec3 SphACS_solar_position(const double T);
	// Calculating the solar position in orthogonal absolute coordinate system
	general::math::Vec3 ACS_solar_position(const double T);

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
	template<class SunlightModel, class ... ArgsType>
	class ISunlight
	{
	public:
		general::math::Vec3 acceleration(
			const general::math::Vec3& solarpos,
			const general::math::Vec3& pointpos,
			const ArgsType ... args) const
		{
			return static_cast<const SunlightModel>(this)->acceleration(solarpos, pointpos, args ...);
		}
	};
}