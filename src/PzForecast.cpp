#include "forecasts/PzForecast.h"
#include "Conversions.h"

namespace ball
{
	PV PzForecast::function(const PV& vec, const JD& t)
	{
		PV func(vec.V1, vec.V2, vec.V3, 0, 0, 0);
		auto xyzPos{ XYZ(vec.P1, vec.P2, vec.P3) };
		double r = xyzPos.length();
		double r_2{ r * r };
		double w_2{ _eW * _eW };
		double h = GCS_height_from_position(xyzPos, _eR, _eFl);
		if (h < MinHeight || h > MaxHeight)
			throw std::runtime_error("Height is out of bounds!");
		// geopotential aceleration with a centrifugal and a coriolis force
		auto xyzAcPot{ _geopotential.acceleration(xyzPos) };
		// atmosphere aceleration a = v * s * rho, 
		// s - a ballistic coefficient,
		// v - a velocity of the vehicle,
		// rho - a density of the atmosphere
		double density = _pAtmosphere->density(xyzPos, t);
		double v = std::sqrt(vec.V1 * vec.V1 + vec.V2 * vec.V2 + vec.V3 * vec.V3);
		double acAtm = v * density * sBall;
		// the addition all the components
		func.V1 = xyzAcPot.X + w_2 * vec.P1 + 2 * _eW * vec.V2 - acAtm * vec.V1;
		func.V2 = xyzAcPot.Y + w_2 * vec.P2 - 2 * _eW * vec.V1 - acAtm * vec.V2;
		func.V3 = xyzAcPot.Z - acAtm * vec.V3;
		return func;
	}
}