#include "RBL.h"

namespace ball
{
	namespace geometry
	{
		RBL& RBL::operator += (const RBL& v) 
		{
			R += v.R;
			B += v.B;
			L += v.L;
			return *this;
		}
		RBL& RBL::operator -= (const RBL& v)
		{
			R -= v.R;
			B -= v.B;
			L -= v.L;
			return *this;
		}
		RBL& RBL::operator /= (const double n)
		{
			R /= n;
			B /= n;
			L /= n;
			return *this;
		}
		RBL& RBL::operator *= (const double n)
		{
			R *= n;
			B *= n;
			L *= n;
			return *this;
		}

		RBL operator + (const RBL& f, const RBL& s)
		{
			return RBL(f.R + s.R, f.B + s.B, f.L + s.L);
		}
		RBL operator - (const RBL& f, const RBL& s)
		{
			return RBL(f.R - s.R, f.B - s.B, f.L - s.L);
		}
		RBL operator * (const RBL& v, const double n)
		{
			return RBL(v.R * n, v.B * n, v.L * n);
		}
		RBL operator / (const RBL& v, const double n)
		{
			return RBL(v.R / n, v.B / n, v.L / n);
		}
		RBL operator * (const double n, const RBL& v)
		{
			return RBL(v.R * n, v.B * n, v.L * n);
		}

		std::ostream& operator << (std::ostream& o, const RBL& v)
		{
			o << "{ " << v.R << "; " << v.B << "; " << v.L << " }";
			return o;
		}
	}
}