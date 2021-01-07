#pragma once
#include <ostream>

namespace ball
{
	namespace geometry
	{
		struct XYZ
		{
			double X, Y, Z;
			XYZ() : X(0), Y(0), Z(0) {}
			XYZ(
				const double x,
				const double y,
				const double z) : X(x), Y(y), Z(z)
			{}
			~XYZ() {}

			double Length() const;

			XYZ& operator += (const XYZ& v);
			XYZ& operator -= (const XYZ& v);
			double operator *= (const XYZ& v);
			XYZ& operator /= (const double n);
			XYZ& operator *= (const double n);

			friend std::ostream& operator << (std::ostream& o, const XYZ& v);

			friend XYZ operator + (const XYZ& f, const XYZ& s);
			friend XYZ operator - (const XYZ& f, const XYZ& s);
			friend XYZ operator * (const XYZ& v, const double n);
			friend XYZ operator / (const XYZ& v, const double n);
			friend XYZ operator * (const double n, const XYZ& v);
		};
	}
}