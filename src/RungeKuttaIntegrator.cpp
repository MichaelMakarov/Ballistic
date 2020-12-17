#include "RungeKuttaIntegrator.h"

namespace ball
{
	types::PV RungeKuttaIntegrator::Integrate(const double h) const
	{
		const double 	h_2 = 0.5 * h,
						h_6 = h / 6;
		types::PV k1 = _func(_x0, _t0),
			k2 = _func(
				types::PV(
					_x0.P1 + h_2 * k1.P1,
					_x0.P2 + h_2 * k1.P2,
					_x0.P3 + h_2 * k1.P3,
					_x0.V1 + h_2 * k1.V1,
					_x0.V2 + h_2 * k1.V2,
					_x0.V3 + h_2 * k1.V3),
				_t0 + h_2),
			k3 = _func(
				types::PV(
					_x0.P1 + h_2 * k2.P1,
					_x0.P2 + h_2 * k2.P2,
					_x0.P3 + h_2 * k2.P3,
					_x0.V1 + h_2 * k2.V1,
					_x0.V2 + h_2 * k2.V2,
					_x0.V3 + h_2 * k2.V3),
				_t0 + h_2),
			k4 = _func(
				types::PV(
					_x0.P1 + h * k3.P1,
					_x0.P2 + h * k3.P2,
					_x0.P3 + h * k3.P3,
					_x0.V1 + h * k3.V1,
					_x0.V2 + h * k3.V2,
					_x0.V3 + h * k3.V3),
				_t0 + h);
		return types::PV(
			_x0.P1 + h_6 * (k1.P1 + 2 * (k2.P1 + k3.P1) + k4.P1),
			_x0.P2 + h_6 * (k1.P2 + 2 * (k2.P2 + k3.P2) + k4.P2),
			_x0.P3 + h_6 * (k1.P3 + 2 * (k2.P3 + k3.P3) + k4.P3),
			_x0.V1 + h_6 * (k1.V1 + 2 * (k2.V1 + k3.V1) + k4.V1),
			_x0.V2 + h_6 * (k1.V2 + 2 * (k2.V2 + k3.V2) + k4.V2),
			_x0.V3 + h_6 * (k1.V3 + 2 * (k2.V3 + k3.V3) + k4.V3)
		);
	}
}