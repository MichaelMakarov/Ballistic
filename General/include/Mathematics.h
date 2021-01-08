#include "GeneralConstants.h"
#include <vector>
#include <ostream>
#include <istream>

namespace ball
{
	namespace math
	{
        inline double DegToRad(const double degrees)
        {
            return degrees * PI / 180.0;
        }
        inline double RadToDeg(const double radians)
        {
            return radians * 180.0 / PI;
        }

        long double Factorial(const size_t x);

		class LegendrePolynomial
		{
		protected:
			size_t _degree;
			std::vector<double> _values;

		private:
			void CreatePolynom(
				std::vector<double>& values,
				const size_t degree);

			bool _normalized;

		public:
			explicit LegendrePolynomial(
				const size_t degree = 0,
				const bool normalized = false);
			~LegendrePolynomial() {}

			size_t Degree() const { return _degree; }

			bool Normalized() const { return _normalized; }

			virtual double operator () (const double x) const;
			double operator [] (const size_t i) const;
		};

		class LegendreFunction : public LegendrePolynomial
		{
		private:
			size_t _derivation;
			bool _normalized;

		public:
			LegendreFunction(
				const size_t degree = 0,
				const size_t derivation = 0,
				const bool normalized = false);
			~LegendreFunction() {}

			size_t Derivation() const { return _derivation; }

			bool Normalized() const { return _normalized; }

			double operator () (const double x) const;
		};

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