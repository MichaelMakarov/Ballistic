#include "LegendrePolynom.h"
#include <cmath>
#include <cstring>
#include "MathFunctions.h"

namespace ball
{
	namespace math
	{
		void LegendrePolynom::CreatePolynom(
			std::vector<double>& values,
			const size_t degree,
			const size_t size)
		{
			values.resize(size);
			std::memset(values.data(), 0, size * sizeof(double));
			switch (degree)
			{
			case 0:
				values[0] = 1.0;
				break;
			case 1:
				values[1] = 1.0;
				break;
			case 2:
				values[0] = -0.5;
				values[2] = 1.5;
				break;
			case 3:
				values[1] = -1.5;
				values[3] = 2.5;
				break;
			case 4:
				values[0] = 3.0 / 8;
				values[2] = -30.0 / 8;
				values[4] = -35.0 / 8;
				break;
			case 5:
				values[1] = 15.0 / 8;
				values[3] = -70.0 / 8;
				values[5] = 63.0 / 8;
				break;
			case 6:
				values[0] = -5.0 / 16;
				values[2] = 105.0 / 16;
				values[4] = -315.0 / 16;
				values[6] = 231.0 / 16;
				break;
			case 7:
				values[1] = -35.0 / 16;
				values[3] = 315.0 / 16;
				values[5] = -693.0 / 16;
				values[7] = 429.0 / 16;
				break;
			case 8:
				values[0] = 35.0 / 128;
				values[2] = -1260.0 / 128;
				values[4] = 6930.0 / 128;
				values[6] = 12012.0 / 128;
				values[8] = 6435.0 / 128;
				break;
			case 9:
				values[1] = 315.0 / 128;
				values[3] = -4620.0 / 128;
				values[5] = 18018.0 / 128;
				values[7] = -25740.0 / 128;
				values[9] = 12155.0 / 128;
				break;
			case 10:
				values[0] = -63.0 / 256;
				values[2] = 3465.0 / 256;
				values[4] = -30030.0 / 256;
				values[6] = 90090.0 / 256;
				values[8] = -109395.0 / 256;
				values[10] = 45189.0 / 256;
				break;
			case 11:
				values[1] = -693.0 / 256;
				values[3] = 15015.0 / 256;
				values[5] = -90090.0 / 256;
				values[7] = 218790.0 / 256;
				values[9] = -230945.0 / 256;
				values[11] = 88179.0 / 256;
				break;
			case 12:
				values[0] = 231.0 / 1024;
				values[2] = -18018.0 / 1024;
				values[4] = 225225.0 / 1024;
				values[6] = -1021020.0 / 1024;
				values[8] = 2078505.0 / 1024;
				values[10] = -1939938.0 / 1024;
				values[12] = 676039.0 / 1024;
				break;
			}
		}

		LegendrePolynom::LegendrePolynom(const size_t degree)
		{
			_norm = std::sqrt(degree + 0.5);
			_degree = degree;
			if (degree <= 12) CreatePolynom(_values, degree, degree + 1);
			else {
				std::vector<double> arr1, arr2;
				CreatePolynom(arr1, 12, degree + 1);
				CreatePolynom(arr2, 11, degree + 1);
				_values.resize(degree + 1);
				for (int n = 12; n < degree; ++n)
				{
					for (int i = n; i >= 0; --i)
						_values[i + 1] = ((2 * n + 1) * arr1[i] - n * arr2[i + 1]) / (n + 1);
					_values[0] = -n * arr2[0] / (n + 1);
					std::memcpy(arr2.data(), arr1.data(), sizeof(double) * (n + 1));
					std::memcpy(arr1.data(), _values.data(), sizeof(double) * (n + 2));
				}
			}
		}

		size_t LegendrePolynom::Degree() const
		{
			return _degree;
		}

		double LegendrePolynom::operator () (const double x) const
		{
			double m = x, s = _values[0];
			for (size_t i = 1; i < _values.size(); ++i)
			{
				s += m * _values[i];
				m *= x;
			}
			return s;
		}

		double LegendrePolynom::operator [] (const size_t i) const
		{
			return _values[i];
		}

		double LegendrePolynom::Norm() const
		{
			return _norm;
		}

		LegendreFunction::LegendreFunction(
			const size_t degree,
			const size_t derivation) : LegendrePolynom(degree)
		{
			_derivation = derivation;
			if (degree < derivation)
			{
				_values.resize(1);
				_values[0] = 0.0;
				_norm = 0.0;
				return;
			}
			else if (derivation > 0)
			{
				for (size_t i = 0; i <= degree - derivation; ++i)
				{
					_values[i] = _values[i + derivation];
					for (size_t m = derivation; m >= 1; --m)
						_values[i] *= i + m;
				}
				_values.resize(degree - derivation + 1);
			}
			long double fact = 1.0;
			for (size_t i = degree - derivation + 1; i <= degree + derivation; ++i)
				fact *= i;
			_norm = std::sqrt((2 * degree + 1) * (derivation == 0 ? 1.0 : 2.0) / fact);
		}

		double LegendreFunction::operator () (const double x) const
		{
			double mult = _derivation > 0 ? std::pow(1 - x * x, _derivation * 0.5) : 1.0;
			return LegendrePolynom::operator() (x) * mult;
		}

		double LegendreFunction::Norm() const
		{
			return _norm;
		}
	}
}