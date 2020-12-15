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
	}
}