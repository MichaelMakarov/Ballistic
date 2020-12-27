#pragma once
#include "DateTime.h"
#include "XYZ.h"

namespace ball
{
    namespace types
    {
        struct StateParams
        {
            geometry::XYZ R;    // a position vector
            geometry::XYZ V;    // a velocity vector
            double Sb;          // a ballistic coefficient
            time::JD T;         // a julian date
            unsigned int LoopN; // a number of the loop

            StateParams(
                const geometry::XYZ& r,
                const geometry::XYZ& v,
                const double sb,
                const time::JD& t,
                const unsigned int n);

            friend std::ostream& operator << (std::ostream& o, const StateParams& p);
        };

        struct OsculParams
        {
            double A;           // a semimajor axis
            double E;           // an eccentricity
            double I;           // an inclination
            double W;           // an argument of periapsis
            double O;           // a longitude of the ascending node
            double M;           // a mean anomaly
            double Sb;          // a ballistic coefficient
            double T;           // a julian date
            unsigned int N;     // a number of the loop

            OsculParams() {}
            OsculParams(
                const double a,
                const double e,
                const double i,
                const double w,
                const double o,
                const double m,
                const double sb,
                const double t,
                const unsigned int vitn = 0);
            ~OsculParams() {}
        };
    }
}