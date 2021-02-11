#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
	namespace space
	{
		// Calculating the solar position in spherical absolute coordinate system
		general::math::Vec3 SphACS_solar_position(const double T);
		// Calculating the solar position in orthogonal absolute coordinate system
		general::math::Vec3 ACS_solar_position(const double T);

		// Solar gravity model
		template<class SolarGravityType>
		class ISolarGravity
		{
		public:
			general::math::Vec3 acceleration(
				const general::math::Vec3& solarpos,
				const general::math::Vec3& pointpos) const
			{
				return static_cast< const SolarGravityType>(this)->acceleration(solarpos, pointpos);
			}
		};

		// Solar pressure model
		template<class SolarLightType, class ... ArgsType>
		class ISolarLight
		{
		public:
			general::math::Vec3 acceleration(
				const general::math::Vec3& solarpos, 
				const general::math::Vec3& pointpos, 
				const ArgsType ... args) const
			{
				return static_cast<const SolarLightType>(this)->acceleration(solarpos, pointpos, args ...);
			}
		};
	}
}