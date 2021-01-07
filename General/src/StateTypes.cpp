#include "StateTypes.h"

namespace ball
{
    namespace types
    {        
        std::ostream& operator << (std::ostream& o, const StateParams& p)
        {
            o << "T: " << p.T.ToDateTime() << "; vec: { " << p.Vec << " }; s = " << p.Sb;
            return o;
        }

        StateParams::StateParams(
            const PV vec,
            const double sb,
            const time::JD& t,
            const size_t loop)
        {
            Vec = vec;
            Sb = sb;
            T = t;
            LoopN = loop;
        }

        StateParams::StateParams(
            const double posX,
            const double posY,
            const double posZ,
            const double velX,
            const double velY,
            const double velZ,
            const double sb,
            const time::JD& t,
            const size_t loop)
        {
            Vec = PV(posX, posY, posZ, velX, velY, velZ);
            Sb = sb;
            T = t;
            LoopN = loop;
        }

        StateParams::StateParams(
            const geometry::XYZ& pos,
            const geometry::XYZ& vel,
            const double sb,
            const time::JD& t,
            const size_t loop)
        {
            Vec = PV(pos, vel);
            Sb = sb;
            T = t;
            LoopN = loop;
        }

        StateParams::StateParams(
            const geometry::RBL& pos,
            const geometry::RBL& vel,
            const double sb,
            const time::JD& t,
            const size_t loop)
        {
            Vec = PV(pos, vel);
            Sb = sb;
            T = t;
            LoopN = loop;
        }

        OsculParams::OsculParams(
            const double a,
            const double e,
            const double i,
            const double w,
            const double o,
            const double m,
            const double sb,
            const time::JD& t,
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
            LoopN = vitn;
        }
    }
}