#pragma once
#include "DateTime.h"
#include "PV.h"

namespace ball
{
    namespace types
    {
        struct StateParams
        {
            PV Vec;             // a vector's coordinates 
            double Sb;          // a ballistic coefficient
            time::JD T;         // a julian date
            size_t LoopN;       // a number of the loop

            StateParams() : Vec(), Sb(0), T(), LoopN(0) {}

            StateParams(
                const PV vec,
                const double sb,
                const time::JD& t,
                const size_t loop);

            StateParams(
                const double posX,
                const double posY,
                const double posZ,
                const double velX,
                const double velY,
                const double velZ,
                const double sb,
                const time::JD& t,
                const size_t loop);

            StateParams(
                const geometry::XYZ& pos,
                const geometry::XYZ& vel,
                const double sb,
                const time::JD& t,
                const size_t loop);

            StateParams(
                const geometry::RBL& pos,
                const geometry::RBL& vel,
                const double sb,
                const time::JD& t,
                const size_t loop);

            friend std::ostream& operator << (std::ostream& o, const StateParams& p);
        };

        struct OsculParams
        {
            double A;               // a semimajor axis
            double E;               // an eccentricity
            double I;               // an inclination
            double W;               // an argument of periapsis
            double O;               // a longitude of the ascending node
            double M;               // a mean anomaly
            double Sb;              // a ballistic coefficient
            time::JD T;               // a julian date
            unsigned int LoopN;     // a number of the loop

            OsculParams() : A(0), E(0), I(0), W(0), O(0), M(0), Sb(0), T(0), LoopN(0) {}
            OsculParams(
                const double a,
                const double e,
                const double i,
                const double w,
                const double o,
                const double m,
                const double sb,
                const time::JD& t,
                const unsigned int vitn = 0);
            ~OsculParams() {}
        };
    }
}