#pragma once
#include "Integrators.h"
#include "general/Matrix.h"

namespace ball
{
	template<size_t n> void fill_array(double(&arr)[n])
	{
		switch (n) {
		case 2: arr[0] = 0.5000; arr[1] = 1; break;
		case 3: arr[0] = 0.2764; arr[1] = 0.7236; arr[2] = 1; break;
		case 4: arr[0] = 0.1727; arr[1] = 0.5000; arr[2] = 0.8273; arr[3] = 1; break;
		case 5: arr[0] = 0.1175; arr[1] = 0.3574; arr[2] = 0.6426; arr[3] = 0.8825; arr[4] = 1; break;
		case 6: arr[0] = 0.0849; arr[1] = 0.2656; arr[2] = 0.5000; arr[3] = 0.7344; arr[4] = 0.9151; arr[5] = 1; break;
		case 7: arr[0] = 0.0641; arr[1] = 0.2041; arr[2] = 0.3954; arr[3] = 0.6046; arr[4] = 0.7959; arr[5] = 0.9359; arr[6] = 1; break;
		}
	}
	template<Arithmetic R, Time T, size_t degree = 4>
	class everhart_integrator : public singlestep_integrator<everhart_integrator<R, T, degree>, R, T>
	{
		static_assert(degree >= 2 && degree <= 7, "Not implemented for Everhart's integrator!");
		
		double _tau[degree]{ 0 };
		general::math::MatrixMxN<degree, degree> _matrix;

	public:
		everhart_integrator() : singlestep_integrator<everhart_integrator<R, T, degree>, R, T>() {
			fill_array(_tau);
			for (size_t i = 0; i < degree; ++i) _matrix(i, i) = 1.0;
			for (size_t i = 1; i < degree; ++i) _matrix(0, i) = -_matrix(0, i - 1) * _tau[i - 1];
			for (size_t i = 1; i < degree; ++i) {
				for (size_t j = i + 1; j < degree; ++j) {
					_matrix(i, j) = _matrix(i - 1, j - 1) - _tau[j - 1] * _matrix(i, j - 1);
				}
			}
		}
		~everhart_integrator() = default;

		template<class Inv> void integrate(
			const R& x0, const T& t0,
			const double step,
			R& xk, T& tk,
			const invoker<Inv, R, const R&, const T&>& func) const
		{
			xk = func(x0, t0);
			auto calc_x = [xk, x0, step](const std::array<R, degree>& c, const double t) {
				R result;
				for (long long i = degree - 1; i >= 0; --i) {
					result += c[i] * (1.0 / (i + 2));
					result *= t;
				}
				return (result + xk) * t * step + x0;
			};
			std::array<R, degree> a, b;
			for (size_t k = 0; k < 3; ++k) {
				// iteration started
				for (size_t i = 0; i < degree; ++i) {
					a[i] = (func(calc_x(b, _tau[i]), t0 + (step * _tau[i])) - xk) / _tau[i];
					for (size_t j = 0; j < i; ++j) {
						a[i] -= a[j];
						a[i] /= (_tau[i] - _tau[j]);
					}
					// recalculating coeffcients of decomposition
					for (size_t m = 0; m < degree; ++m) {
						b[m] = R();
						for (size_t n = m; n < degree; ++n)
							b[m] += a[n] * _matrix(m, n);
					}
				}
			}
			tk = t0 + step;
			for (size_t i = 0; i < degree; ++i)
				xk += b[i] * (1.0 / (i + 2));
			xk *= step;
			xk += x0;
		}
	};

}