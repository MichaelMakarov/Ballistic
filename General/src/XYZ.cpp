#include "XYZ.h"
#include <cmath>

namespace ball
{
	namespace geometry
	{
		double XYZ::Length() const
		{
			return std::sqrt(X * X + Y * Y + Z * Z);
		}

        XYZ operator + (const XYZ& f, const XYZ& s)
        {
            return XYZ(f.X + s.X, f.Y + s.Y, f.Z + s.Z);
        }
        XYZ operator - (const XYZ& f, const XYZ& s)
        {
            return XYZ(f.X - s.X, f.Y - s.Y, f.Z - s.Z);
        }
        XYZ operator * (const XYZ& v, const double n)
        {
            return XYZ(v.X * n, v.Y * n, v.Z * n);
        }
        XYZ operator / (const XYZ& v, const double n)
        {
            return XYZ(v.X / n, v.Y / n, v.Z / n);
        }
        XYZ operator * (const double n, const XYZ& v)
        {
            return XYZ(v.X * n, v.Y * n, v.Z * n);
        }

        std::ostream& operator << (std::ostream& o, const XYZ& v)
        {
            o << "{ " << v.X << "; " << v.Y << "; " << v.Z << " }";
            return o;
        }
	}
}