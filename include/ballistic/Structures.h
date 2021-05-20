#pragma once
#include "general/Vector.h"
#include "general/Times.h"
#include <istream>

namespace ball
{
    /// <summary>
    /// 6d vector (x, y, z, vx, vy, vz)
    /// </summary>
    using Vec6 = general::math::Vec<6>;

    // Struct represents the movement of the center of mass.
    // Contains major parameters of the movement.
    struct State
    {
        Vec6 vec;                   // vector's coordinates (depends on coordinate system)
        double Sb;                  // a ballistic coefficient
        general::time::JD T;        // a julian date
        size_t loop;                // a number of the loop

        State() : vec(), Sb(0), T(), loop(0) {}
        State(
            const Vec6 vec,
            const double sb,
            const general::time::JD& t,
            const size_t loop) :
            vec{ vec }, Sb{ sb }, T{ t }, loop{ loop }
        {}
        State(
            const double posX,
            const double posY,
            const double posZ,
            const double velX,
            const double velY,
            const double velZ,
            const double sb,
            const general::time::JD& t,
            const size_t loop) :
            vec{ { posX, posY, posZ, velX, velY, velZ } },
            Sb{ sb }, T{ t }, loop{ loop }
        {}

        friend std::ostream& operator << (std::ostream& os, const State& p);
    };

    // Struct represents the oscullar parameters of the orbit.
    // Contains the oscullar parameters, time, ballstic coefficient and loop.
    struct Oscul
    {
        double A;               // a semimajor axis
        double E;               // an eccentricity
        double I;               // an inclination
        double W;               // an argument of periapsis
        double O;               // a longitude of the ascending node
        double M;               // a mean anomaly
        double Sb;              // a ballistic coefficient
        general::time::JD T;    // a julian date refered to midnight
        size_t loop;            // a number of the loop

        Oscul() : A(0), E(0), I(0), W(0), O(0), M(0), Sb(0), T(0), loop(0) {}
        Oscul(
            const double a,
            const double e,
            const double i,
            const double w,
            const double o,
            const double m,
            const double sb,
            const general::time::JD& t,
            const size_t loop = 0) :
            A{ a }, E{ e }, I{ i }, W{ w }, O{ o }, M{ m }, Sb{ sb }, T{ t }, loop{ loop }
        {}
        ~Oscul() {}

        friend std::ostream& operator << (std::ostream& os, const Oscul& p);
    };
}