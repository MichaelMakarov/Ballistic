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

			friend std::ostream& operator << (std::ostream& o, const RBL& v);

		};
	}
}