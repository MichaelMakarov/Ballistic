#pragma once
#include "Forecast.h"

namespace ball
{
	// A forecast implements solar radiation pressure
	class AngForecast : public Forecast<AngForecast, general::math::PV>
	{
	private:


	public:
		// Calculating accelerations
		general::math::PV function(
			const general::math::PV& rot, 
			const general::time::JD& t, 
			const general::math::PV& pos);
	};
}