#include "PlahovIntegrator.h"
#include <algorithm>

namespace ball
{
	types::PV PlahovIntegrator::Integrate(const double step) const
	{
		const size_t n = 3;
		double h = step,
			d = 10.0,
			eps = 1e-3;
		types::PV fx0 = _func(_x0, _t0), fx[3];
		math::Matrix3x3 pqr, ref = math::Matrix3x3::Eye();
		types::PV xi;
		while (d > eps)
		{
			ref = pqr;
			for (size_t i = 0; i < n; ++i)
			{
				xi = types::PV(
					_x0.P1 + _a[i] * h * (_x0.V1 + _a[i] * h * (0.5 * fx0.V1 + _a[i] * (1.0 / 6 * pqr(0, 0) + _a[i] * (1.0 / 12 * pqr(1, 0) + 0.05 * _a[i] * pqr(2, 0))))),
					_x0.P2 + _a[i] * h * (_x0.V2 + _a[i] * h * (0.5 * fx0.V2 + _a[i] * (1.0 / 6 * pqr(0, 1) + _a[i] * (1.0 / 12 * pqr(1, 1) + 0.05 * _a[i] * pqr(2, 1))))),
					_x0.P3 + _a[i] * h * (_x0.V3 + _a[i] * h * (0.5 * fx0.V3 + _a[i] * (1.0 / 6 * pqr(0, 2) + _a[i] * (1.0 / 12 * pqr(1, 2) + 0.05 * _a[i] * pqr(2, 2))))),
					_x0.V1 + _a[i] * h * (fx0.V1 + _a[i] * (0.5 * pqr(0, 0) + _a[i] * (1.0 / 3 * pqr(1, 0) + 0.25 * _a[i] * pqr(2, 0)))),
					_x0.V2 + _a[i] * h * (fx0.V2 + _a[i] * (0.5 * pqr(0, 1) + _a[i] * (1.0 / 3 * pqr(1, 1) + 0.25 * _a[i] * pqr(2, 1)))),
					_x0.V3 + _a[i] * h * (fx0.V3 + _a[i] * (0.5 * pqr(0, 2) + _a[i] * (1.0 / 3 * pqr(1, 2) + 0.25 * _a[i] * pqr(2, 2))))
				);
				fx[i] = _func(xi, _t0 + _a[i] * h);
			}
			pqr = _mA * math::Matrix3x3(
				fx[0].V1 - fx0.V1, fx[0].V2 - fx0.V2, fx[0].V3 - fx0.V3,
				fx[1].V2 - fx0.V1, fx[1].V2 - fx0.V2, fx[1].V3 - fx0.V3,
				fx[2].V3 - fx0.V1, fx[2].V2 - fx0.V2, fx[2].V3 - fx0.V3);
			d = std::max({ 
				std::abs(pqr[0] - ref[0]), 
				std::abs(pqr[1] - ref[1]),
				std::abs(pqr[2] - ref[2]),
				std::abs(pqr[3] - ref[3]),
				std::abs(pqr[4] - ref[4]),
				std::abs(pqr[5] - ref[5]),
				std::abs(pqr[6] - ref[6]),
				std::abs(pqr[7] - ref[7]),
				std::abs(pqr[8] - ref[8])
			});
		}
		return types::PV(
			_x0.P1 + h * (_x0.V1 + h * (0.5 * fx0.V1 + 1.0 / 6 * pqr(0, 0) + 1.0 / 12 * pqr(1, 0) + 0.05 * pqr(2, 0))),
			_x0.P2 + h * (_x0.V2 + h * (0.5 * fx0.V2 + 1.0 / 6 * pqr(0, 1) + 1.0 / 12 * pqr(1, 1) + 0.05 * pqr(2, 1))),
			_x0.P3 + h * (_x0.V3 + h * (0.5 * fx0.V3 + 1.0 / 6 * pqr(0, 2) + 1.0 / 12 * pqr(1, 2) + 0.05 * pqr(2, 2))),
			_x0.V1 + h * (fx0.V1 + 0.5 * pqr(0, 0) + 1.0 / 3 * pqr(1, 0) + 0.25 * pqr(2, 0)),
			_x0.V2 + h * (fx0.V2 + 0.5 * pqr(0, 1) + 1.0 / 3 * pqr(1, 1) + 0.25 * pqr(2, 1)),
			_x0.V3 + h * (fx0.V3 + 0.5 * pqr(0, 2) + 1.0 / 3 * pqr(1, 2) + 0.25 * pqr(2, 2))
		);
	}
}