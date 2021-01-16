#include "Geometry.h"
#include <cmath>

namespace ball
{
    namespace geometry
    {
        double XYZ::length() const
        {
            return std::sqrt(X * X + Y * Y + Z * Z);
        }

        XYZ& XYZ::operator += (const XYZ& v)
        {
            X += v.X;
            Y += v.Y;
            Z += v.Z;
            return *this;
        }
        XYZ& XYZ::operator -= (const XYZ& v)
        {
            X -= v.X;
            Y -= v.Y;
            Z -= v.Z;
            return *this;
        }
        double XYZ::operator *= (const XYZ& v)
        {
            return X * v.X + Y * v.Y + Z * v.Z;
        }
        XYZ& XYZ::operator *= (const double n)
        {
            X *= n;
            Y *= n;
            Z *= n;
            return *this;
        }
        XYZ& XYZ::operator /= (const double n)
        {
            X /= n;
            Y /= n;
            Z /= n;
            return *this;
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

        std::ostream& operator << (std::ostream& o, const PV& pv)
        {
            o << pv.P1 << "; " << pv.P2 << "; " << pv.P3 << "; " <<
                pv.V1 << "; " << pv.V2 << "; " << pv.V3;
            return o;
        }

        PV& PV::operator += (const PV& pv)
        {
            P1 += pv.P1;
            P2 += pv.P2;
            P3 += pv.P3;
            V1 += pv.V1;
            V2 += pv.V2;
            V3 += pv.V3;
            return *this;
        }
        PV& PV::operator -= (const PV& pv)
        {
            P1 -= pv.P1;
            P2 -= pv.P2;
            P3 -= pv.P3;
            V1 -= pv.V1;
            V2 -= pv.V2;
            V3 -= pv.V3;
            return *this;
        }
        PV& PV::operator *= (const double m)
        {
            P1 *= m;
            P2 *= m;
            P3 *= m;
            V1 *= m;
            V2 *= m;
            V3 *= m;
            return *this;
        }

        PV operator + (const PV& f, const PV& s)
        {
            return PV(
                f.P1 + s.P1,
                f.P2 + s.P2,
                f.P3 + s.P3,
                f.V1 + s.V1,
                f.V2 + s.V2,
                f.V3 + s.V3
            );
        }
        PV operator - (const PV& f, const PV& s)
        {
            return PV(
                f.P1 - s.P1,
                f.P2 - s.P2,
                f.P3 - s.P3,
                f.V1 - s.V1,
                f.V2 - s.V2,
                f.V3 - s.V3
            );
        }
        PV operator * (const PV& pv, const double m)
        {
            return PV(
                pv.P1 * m,
                pv.P2 * m,
                pv.P3 * m,
                pv.V1 * m,
                pv.V2 * m,
                pv.V3 * m
            );
        }
        PV operator * (const double m, const PV& pv)
        {
            return PV(
                pv.P1 * m,
                pv.P2 * m,
                pv.P3 * m,
                pv.V1 * m,
                pv.V2 * m,
                pv.V3 * m
            );
        }
    }
}