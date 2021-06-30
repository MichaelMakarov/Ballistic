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

    // The parameters of motion
    // Contains major parameters of the movement.
    template<size_t dim = 6>
    struct Params
    {
        general::math::Vec<dim> vec{};// vector defines an orbit movement
        general::time::JD T{};        // a julian date
        double sb{};                  // ballistic parameter
        size_t loop{};                // a number of the loop

        friend std::ostream& operator << (std::ostream& os, const Params& p)
        {
            os << "T: " << jd_to_datetime(p.T) << "; vec: " << p.vec << "; loop = " << p.loop << "; s = " << p.sb;
            return os;
        }
    };

    // Struct represents the oscullar parameters of the orbit.
    // Contains the oscullar parameters, time, ballstic coefficient and loop.
    struct Oscul
    {
        double semiaxis;                // a semimajor axis
        double eccentricity;            // an eccentricity
        double inclination;             // an inclination
        double periapsis;               // an argument of periapsis
        double ascendnode;              // a longitude of the ascending node
        double meananomaly;             // a mean anomaly
        double ecanomaly;               // an eccentric anomaly
        double trueanomaly;             // a true anomaly
        double latitudearg;                 // a longitude argument
    };
}