#pragma once
#include <ostream>

namespace ball
{
	namespace geometry
	{
		struct RBL
		{
			double R, B, L;

			RBL() : R(0), B(0), L(0) {}
			RBL(
				const double r,
				const double b,
				const double l) : R(r), B(b), L(l)
			{}
			~RBL() {}

			RBL& operator += (const RBL& v);
			RBL& operator -= (const RBL& v);
			RBL& operator /= (const double n);
			RBL& operator *= (const double n);

			friend std::ostream& operator << (std::ostream& o, const RBL& v);

			friend RBL operator + (const RBL& f, const RBL& s);
			friend RBL operator - (const RBL& f, const RBL& s);
			friend RBL operator * (const RBL& v, const double n);
			friend RBL operator / (const RBL& v, const double n);
			friend RBL operator * (const double n, const RBL& v);

		};
	}
}