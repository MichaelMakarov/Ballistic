#include "LegendrePolynom.h"
#include <cmath>
#include <cstring>
#include "MathFunctions.h"

namespace ball
{
	namespace math
	{
		void LegendrePolynomial::CreatePolynom(
			std::vector<double>& values,
			const size_t degree)
		{
			size_t size = degree + 1;
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
				values[4] = 35.0 / 8;
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
				values[6] = -12012.0 / 128;
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
				values[10] = 46189.0 / 256;
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
			case 13:
				values[1] = 3003.0 / 1024;
				values[3] = -90090.0 / 1024;
				values[5] = 765765.0 / 1024;
				values[7] = -2771340.0 / 1024;
				values[9] = 4849845.0 / 1024;
				values[11] = -4056234.0 / 1024;
				values[13] = 1300075.0 / 1024;
				break;
			case 14:
				values[0] = -429.0 / 2048;
				values[2] = 45045.0 / 2048;
				values[4] = -765765.0 / 2048;
				values[6] = 4849845.0 / 2048;
				values[8] = -14549535.0 / 2048;
				values[10] = 22309287.0 / 2048;
				values[12] = -16900975.0 / 2048;
				values[14] = 5014575.0 / 2048;
				break;
			case 15:
				values[1] = -6435.0 / 2048;
				values[3] = 255255.0 / 2048;
				values[5] = -2909907.0 / 2048;
				values[7] = 14549535.0 / 2048;
				values[9] = -37182145.0 / 2048;
				values[11] = 50702925.0 / 2048;
				values[13] = -35102025.0 / 2048;
				values[15] = 9694845.0 / 2048;
				break;
			case 16:
				values[0] = 6435.0 / 32768;
				values[2] = -875160.0 / 32768;
				values[4] = 19399380.0 / 32768;
				values[6] = -162954792.0 / 32768;
				values[8] = 669278610.0 / 32768;
				values[10] = -1487285800.0 / 32768;
				values[12] = 1825305300.0 / 32768;
				values[14] = -1163381400.0 / 32768;
				values[16] = 300540195.0 / 32768;
				break;
			case 17:
				values[1] = 109395.0 / 32768;
				values[3] = -5542680.0 / 32768;
				values[5] = 81477396.0 / 32768;
				values[7] = -535422888.0 / 32768;
				values[9] = 1859107250.0 / 32768;
				values[11] = -3650610600.0 / 32768;
				values[13] = 4071834900.0 / 32768;
				values[15] = -2404321560.0 / 32768;
				values[17] = 583401555.0 / 32768;
				break;
			}
		}

		LegendrePolynomial::LegendrePolynomial(
			const size_t degree,
			const bool normalized)
		{
			_normalized = normalized;
			_degree = degree;
			if (degree <= 17) CreatePolynom(_values, degree);
			else {
				std::vector<double> arr1, arr2;
				CreatePolynom(arr1, 17);
				CreatePolynom(arr2, 16);
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
			if (normalized)
			{
				const double norm = std::sqrt(degree + 0.5);
				for (double& v : _values) v *= norm;
			}
		}

		size_t LegendrePolynomial::Degree() const
		{
			return _degree;
		}

		bool LegendrePolynomial::Normalized() const
		{
			return _normalized;
		}

		double LegendrePolynomial::operator () (const double x) const
		{
			double m = x, s = _values[0];
			for (size_t i = 1; i < _values.size(); ++i)
			{
				s += m * _values[i];
				m *= x;
			}
			return s;
		}

		double LegendrePolynomial::operator [] (const size_t i) const
		{
			return _values[i];
		}

		LegendreFunction::LegendreFunction(
			const size_t degree,
			const size_t derivation,
			const bool normalized) : LegendrePolynomial(degree)
		{
			_derivation = derivation;
			_normalized = normalized;
			if (degree < derivation)
			{
				_values.resize(1);
				_values[0] = 0.0;
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
			if (normalized)
			{
				long double fact = 1.0;
				for (size_t i = degree - derivation + 1; i <= degree + derivation; ++i)
					fact *= i;
				double norm = std::sqrt((2 * degree + 1) * (derivation == 0 ? 1.0 : 2.0) / fact);
				for (double& v : _values) v *= norm;
			}
		}

		double LegendreFunction::operator () (const double x) const
		{
			double mult = _derivation > 0 ? std::pow(1 - x * x, _derivation * 0.5) : 1.0;
			return LegendrePolynomial::operator() (x) * mult;
		}

		bool LegendreFunction::Normalized() const
		{
			return _normalized;
		}
	}
}