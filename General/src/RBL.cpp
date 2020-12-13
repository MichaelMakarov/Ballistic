#include "RBL.h"

namespace ball
{
	namespace geometry
	{
		std::ostream& operator << (std::ostream& o, const RBL& v)
		{
			o << "{ " << v.R << "; " << v.B << "; " << v.L << " }";
			return o;
		}
	}
}