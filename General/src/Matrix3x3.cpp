#include "Matrix3x3.h"
#include <memory>

namespace ball
{
	namespace math
	{
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

		Matrix3x3 Matrix3x3::Inv(const Matrix3x3& m)
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
			std::_Xruntime_error("Degenerate matrix!");
		}

		Matrix3x3 Matrix3x3::Eye()
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
	}
}