#pragma once
#include "general/Geometry.h"
#include "general/Times.h"
#include <istream>

namespace ball
{
    // Struct represents the movement of the center of mass.
    // Contains major parameters of the movement.
    struct State
    {
        general::math::PV Vec;      // vector's coordinates (depends on coordinate system)
        double Sb;                  // a ballistic coefficient
        general::time::JD T;        // a julian date
        size_t Loop;                // a number of the loop

        State() : Vec(), Sb(0), T(), Loop(0) {}
        State(
            const general::math::PV vec,
            const double sb,
            const general::time::JD& t,
            const size_t loop) :
            Vec{ vec }, Sb{ sb }, T{ t }, Loop{ loop }
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
            Vec{ posX, posY, posZ, velX, velY, velZ },
            Sb{ sb }, T{ t }, Loop{ loop }
        {}
        State(
            const general::math::Vec3& pos,
            const general::math::Vec3& vel,
            const double sb,
            const general::time::JD& t,
            const size_t loop) :
            Vec{ general::math::PV(pos, vel) }, Sb{ sb }, T{ t }, Loop{ loop }
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
        general::time::JD T;             // a julian date
        size_t Loop;            // a number of the loop

        Oscul() : A(0), E(0), I(0), W(0), O(0), M(0), Sb(0), T(0), Loop(0) {}
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
            A{ a }, E{ e }, I{ i }, W{ w }, O{ o }, M{ m }, Sb{ sb }, T{ t }, Loop{ loop }
        {}
        ~Oscul() {}

        friend std::ostream& operator << (std::ostream& os, const Oscul& p);
    };

    
}