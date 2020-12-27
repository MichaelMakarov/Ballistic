#include "RungeKuttaIntegrator.h"

namespace ball
{
	types::PV RungeKuttaIntegrator::Integrate(const double h) const
	{
		const double 	h_2 = 0.5 * h,
						h_6 = h / 6;
		types::PV k1 = _func(_x0, _t0),
			k2 = _func(_x0 + h_2 * k1, _t0 + h_2),
			k3 = _func(_x0 + h_2 * k2, _t0 + h_2),
			k4 = _func(_x0 + h 	 * k3, _t0 + h);
		return _x0 + h_6 * (k1 + 2 * (k2 + k3) + k4);
	}
}