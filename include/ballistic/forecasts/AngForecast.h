#pragma once
#include "Forecast.h"

namespace ball
{
	// A forecast implements solar radiation pressure
	class AngForecast : public Forecast<AngForecast>
	{
	private:


	public:
		// Calculating accelerations
		general::math::PV function(const general::math::PV& vec, const general::time::JD& t);
	};
}