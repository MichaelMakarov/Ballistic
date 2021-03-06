#pragma once
#include "SolarModel.h"

namespace ball
{
	class SunlightPressure : public ISunlight<SunlightPressure, double>
	{
	private:
		const double _q{ 4.65e5 };

	public:
		general::math::Vec3 acceleration(
			const general::math::Vec3& solarpos,
			const general::math::Vec3& pointpos,
			const double k) const;
	};
}