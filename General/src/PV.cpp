#include "PV.h"

namespace ball
{
	namespace types
	{
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
