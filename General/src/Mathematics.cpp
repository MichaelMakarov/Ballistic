#include "Mathematics.h"

namespace ball
{
	namespace math
	{
		void LegendrePolynomial::create(
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
			if (degree <= 17) create(_values, degree);
			else {
				std::vector<double> arr1, arr2;
				create(arr1, 17);
				create(arr2, 16);
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

		Matrix3x3::Matrix3x3(const Matrix3x3& m)
		{
			std::memcpy(_values, m._values, sizeof(double) * 9);
		}

		Matrix3x3& Matrix3x3::operator = (const Matrix3x3& m)
		{
			std::memcpy(_values, m._values, sizeof(double) * 9);
			return *this;
		}

		double Matrix3x3::Det() const
		{
			return	_values[0] * (_values[4] * _values[8] - _values[5] * _values[7]) -
				_values[1] * (_values[3] * _values[8] - _values[5] * _values[6]) +
				_values[2] * (_values[3] * _values[7] - _values[4] * _values[6]);
		}

		double Matrix3x3::operator () (const size_t m, const size_t n) const
		{
			return _values[m * 3 + n];
		}

		double Matrix3x3::operator [] (const size_t i) const
		{
			return _values[i];
		}

		Matrix3x3 Matrix3x3::inv(const Matrix3x3& m)
		{
			double det = m.Det();
			if (std::abs(det) > 1e-16)
				return Matrix3x3(
					(m._values[4] * m._values[8] - m._values[5] * m._values[7]) / det,
					(m._values[2] * m._values[7] - m._values[1] * m._values[8]) / det,
					(m._values[1] * m._values[5] - m._values[2] * m._values[4]) / det,
					(m._values[5] * m._values[6] - m._values[3] * m._values[8]) / det,
					(m._values[0] * m._values[8] - m._values[2] * m._values[6]) / det,
					(m._values[2] * m._values[3] - m._values[0] * m._values[5]) / det,
					(m._values[3] * m._values[7] - m._values[4] * m._values[6]) / det,
					(m._values[1] * m._values[6] - m._values[0] * m._values[7]) / det,
					(m._values[0] * m._values[4] - m._values[1] * m._values[3]) / det
				);
			throw std::runtime_error("Degenerate matrix!");
		}

		Matrix3x3 Matrix3x3::eye()
		{
			return Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
		}

		Matrix3x3& Matrix3x3::operator += (const Matrix3x3& m)
		{
			for (size_t i = 0; i < 9; ++i)
				_values[i] += m._values[i];
			return *this;
		}
		Matrix3x3& Matrix3x3::operator -= (const Matrix3x3& m)
		{
			for (size_t i = 0; i < 9; ++i)
				_values[i] -= m._values[i];
			return *this;
		}
		Matrix3x3& Matrix3x3::operator *= (const Matrix3x3& m)
		{
			double buf[9]{ 0 };
			for (size_t i = 0; i < 3; ++i)
				for (size_t j = 0; j < 3; ++j)
					for (size_t k = 0; k < 3; ++k)
						buf[i * 3 + j] += _values[i * 3 + k] * m._values[k * 3 + j];
			std::memcpy(_values, buf, sizeof(double) * 9);
			return *this;
		}
		Matrix3x3& Matrix3x3::operator *= (const double v)
		{
			for (size_t i = 0; i < 9; ++i)
				_values[i] *= v;
			return *this;
		}
		Matrix3x3& Matrix3x3::operator /= (const double v)
		{
			for (size_t i = 0; i < 9; ++i)
				_values[i] /= v;
			return *this;
		}

		Matrix3x3 operator + (const Matrix3x3& m1, const Matrix3x3& m2)
		{
			return Matrix3x3(
				m1._values[0] + m2._values[0],
				m1._values[1] + m2._values[1],
				m1._values[2] + m2._values[2],
				m1._values[3] + m2._values[3],
				m1._values[4] + m2._values[4],
				m1._values[5] + m2._values[5],
				m1._values[6] + m2._values[6],
				m1._values[7] + m2._values[7],
				m1._values[8] + m2._values[8]
			);
		}
		Matrix3x3 operator - (const Matrix3x3& m1, const Matrix3x3& m2)
		{
			return Matrix3x3(
				m1._values[0] - m2._values[0],
				m1._values[1] - m2._values[1],
				m1._values[2] - m2._values[2],
				m1._values[3] - m2._values[3],
				m1._values[4] - m2._values[4],
				m1._values[5] - m2._values[5],
				m1._values[6] - m2._values[6],
				m1._values[7] - m2._values[7],
				m1._values[8] - m2._values[8]
			);
		}
		Matrix3x3 operator * (const Matrix3x3& m1, const Matrix3x3& m2)
		{
			double buf[9]{ 0 };
			for (size_t i = 0; i < 3; ++i)
				for (size_t j = 0; j < 3; ++j)
					for (size_t k = 0; k < 3; ++k)
						buf[i * 3 + j] += m1._values[i * 3 + k] * m2._values[k * 3 + j];
			return Matrix3x3(
				buf[0], buf[1], buf[2],
				buf[3], buf[4], buf[5],
				buf[6], buf[7], buf[8]);
		}
		Matrix3x3 operator * (const Matrix3x3& m, const double v)
		{
			return Matrix3x3(
				m._values[0] * v,
				m._values[1] * v,
				m._values[2] * v,
				m._values[3] * v,
				m._values[4] * v,
				m._values[5] * v,
				m._values[6] * v,
				m._values[7] * v,
				m._values[8] * v
			);
		}
		Matrix3x3 operator * (const double v, const Matrix3x3& m)
		{
			return Matrix3x3(
				m._values[0] * v,
				m._values[1] * v,
				m._values[2] * v,
				m._values[3] * v,
				m._values[4] * v,
				m._values[5] * v,
				m._values[6] * v,
				m._values[7] * v,
				m._values[8] * v
			);
		}
		Matrix3x3 operator / (const Matrix3x3& m, const double v)
		{
			return Matrix3x3(
				m._values[0] / v,
				m._values[1] / v,
				m._values[2] / v,
				m._values[3] / v,
				m._values[4] / v,
				m._values[5] / v,
				m._values[6] / v,
				m._values[7] / v,
				m._values[8] / v
			);
		}

		std::ostream& operator << (std::ostream& o, const Matrix3x3& m)
		{
			o << "{ { " << m._values[0] << "; " << m._values[1] << "; " << m._values[2] <<
				" } { " << m._values[3] << "; " << m._values[4] << "; " << m._values[5] <<
				" } { " << m._values[6] << "; " << m._values[7] << "; " << m._values[8] <<
				" } }";
			return o;
		}
		std::istream& operator >> (std::istream& is, Matrix3x3& m)
		{
			for (size_t i = 0; i < 9; ++i)
				is >> m._values[i];
			return is;
		}

        long double factorial(const size_t x)
        {
            static const long double values[41]
            {
                1,
                1,
                2,
                6,
                24,
                120,
                720,
                5040,
                40320,
                362880,
                3628800,
                39916800,
                479001600,
                6227020800,
                87178291200,
                1307674368000,
                20922789888000,
                355687428096000,
                6402373705728000,
                121645100408832000,
                2432902008176640000.0,
                5.109094217170944e+19L,
                1.1240007277776077e+21L,
                2.5852016738884978e+22L,
                6.2044840173323941e+23L,
                1.5511210043330986e+25L,
                4.0329146112660565e+26L,
                1.0888869450418352e+28L,
                3.0488834461171384e+29L,
                8.8417619937397008e+30L,
                2.6525285981219103e+32L,
                8.2228386541779224e+33L,
                2.6313083693369352e+35L,
                8.6833176188118859e+36L,
                2.9523279903960412e+38L,
                1.0628380765425749e+40L,
                3.8262170755532692e+41L,
                1.4157003179547095e+43L,
                5.3796612082278958e+44L,
                2.0980678712088794e+46L,
                8.3922714848355173e+45L
            };
            if (x < 40) return values[x];
            throw std::runtime_error("Unsupported!");
        }
	}
}