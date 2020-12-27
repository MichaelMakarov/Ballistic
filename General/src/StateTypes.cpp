#include "StateTypes.h"

namespace ball
{
    namespace types
    {        
        std::ostream& operator << (std::ostream& o, const StateParams& p)
        {
            o << "T: " << p.T.ToDateTime() << "r: " << p.R << " v: " << p.V << " s = " << p.Sb;
            return o;
        }

        StateParams::StateParams(
            const geometry::XYZ& r,
            const geometry::XYZ& v,
            const double sb,
            const time::JD& t,
            const unsigned int n)
        {
            R = r;
            V = v;
            Sb = sb;
            T = t;
            LoopN = n;
        }

        OsculParams::OsculParams(
            const double a,
            const double e,
            const double i,
            const double w,
            const double o,
            const double m,
            const double sb,
            const double t,
            const unsigned int vitn)
        {
            A = a;
            E = e;
            I = i;
            W = w;
            O = o;
            M = m;
            Sb = sb;
            T = t;
            N = vitn;
        }
    }
}