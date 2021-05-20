#include "Structures.h"
#include "general/GeneralConstants.h"
#include "Conversions.h"

namespace ball
{
    std::ostream& operator<<(std::ostream& os, const State& p)
    {
        os << "T: " << p.T.to_datetime() << "; Vec: " << p.vec << "; s = " <<
            p.Sb << "; loop = " << p.loop;
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const Oscul& p)
    {
        os << "T: " << p.T << "; a = " << p.A << "; e = " << p.E <<
            "; i = " << p.I << "; w = " << p.W << "; o = " << p.O <<
            "; M = " << p.M << "; s = " << p.Sb << "; loop = " << p.loop;
        return os;
    }

	
}