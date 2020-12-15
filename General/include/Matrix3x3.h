#pragma once
#include <ostream>
#include <istream>

namespace ball 
{
	namespace math
	{
		class Matrix3x3
		{
		private:
			double _values[9];

		public:
			Matrix3x3() : _values{ 0, 0, 0, 0, 0, 0, 0, 0, 0 }
			{}
			Matrix3x3(
				const double m11, const double m12, const double m13,
				const double m21, const double m22, const double m23,
				const double m31, const double m32, const double m33)
				: _values{ m11, m12, m13, m21, m22, m23, m31, m32, m33 }
			{}
			Matrix3x3(const Matrix3x3& m);
			~Matrix3x3() {}

			Matrix3x3& operator = (const Matrix3x3& m);

			double Det() const;

			double operator () (const size_t m, const size_t n) const;
			double operator [] (const size_t i) const;
			
			static Matrix3x3 Inv(const Matrix3x3& m);
			static Matrix3x3 Eye();

			Matrix3x3& operator += (const Matrix3x3& m);
			Matrix3x3& operator -= (const Matrix3x3& m);
			Matrix3x3& operator *= (const Matrix3x3& m);
			Matrix3x3& operator *= (const double v);
			Matrix3x3& operator /= (const double v);

			friend Matrix3x3 operator + (const Matrix3x3& m1, const Matrix3x3& m2);
			friend Matrix3x3 operator - (const Matrix3x3& m1, const Matrix3x3& m2);
			friend Matrix3x3 operator * (const Matrix3x3& m1, const Matrix3x3& m2);
			friend Matrix3x3 operator * (const Matrix3x3& m, const double v);
			friend Matrix3x3 operator / (const Matrix3x3& m, const double v);
			friend Matrix3x3 operator * (const double v, const Matrix3x3& m);

			friend std::ostream& operator << (std::ostream& o, const Matrix3x3& m);
			friend std::istream& operator >> (std::istream& i, Matrix3x3& m);
		};
	}
}