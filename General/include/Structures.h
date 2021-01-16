#pragma once
#include "Geometry.h"
#include "Times.h"
#include <istream>

namespace ball
{
	namespace geometry
	{
        // Struct represents the movement of the center of mass.
        // Contains major parameters of the movement.
		struct State
		{
            PV Vec;             // vector's coordinates (depends on coordinate system)
            double Sb;          // a ballistic coefficient
            time::JD T;         // a julian date
            size_t Loop;        // a number of the loop

            State() : Vec(), Sb(0), T(), Loop(0) {}
            State(
                const PV vec,
                const double sb,
                const time::JD& t,
                const size_t loop) :
                Vec{ vec }, Sb{ sb }, T{ t }, Loop{loop}
            {}
            State(
                const double posX,
                const double posY,
                const double posZ,
                const double velX,
                const double velY,
                const double velZ,
                const double sb,
                const time::JD& t,
                const size_t loop) :
                Vec{ posX, posY, posZ, velX, velY, velZ },
                Sb{ sb }, T{ t }, Loop{ loop }
            {}
            State(
                const geometry::XYZ& pos,
                const geometry::XYZ& vel,
                const double sb,
                const time::JD& t,
                const size_t loop) :
                Vec{ PV(pos, vel) }, Sb{ sb }, T{ t }, Loop{ loop }
            {}
            State(
                const geometry::RBL& pos,
                const geometry::RBL& vel,
                const double sb,
                const time::JD& t,
                const size_t loop) :
                Vec{ PV(pos, vel) }, Sb{ sb }, T{ t }, Loop{ loop } 
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
            time::JD T;             // a julian date
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
                const time::JD& t,
                const size_t loop = 0) :
                A{ a }, E{ e }, I{ i }, W{ w }, O{ o }, M{ m }, Sb{ sb }, T{ t }, Loop{ loop }
            {}
            ~Oscul() {}

            friend std::ostream& operator << (std::ostream& os, const Oscul& p);
        };

        // Struct represents the tle parameters of the movement of the center of mass.
        struct TLE
        {
            std::string Header;             // a header of the file

            unsigned int Satellite;         // a number of the sattelite in the NORAD database
            char Classification;            // classification (U = Unclassified - i.e. not secret)
            unsigned int LaunchYear;        // two last digits of the year when the sattelite was launched
            unsigned int LaunchNumber;      // a launch number of the year
            char LaunchPlace[4] = { 0 };    // a launch part (place)
            unsigned int EpochYear;         // a year of the epoch (two last digits)
            double EpochTime;               // time of the epoch (a day number of the year with a part of the day)
            double FirstAcc;                // a first derivative of the movement (an acceleration, loop per day^2)
            double SecondAcc;               // a second derivative of the movement divided by 6 (loop per day^3)
            double Deceleration;            // a deceleration coefficient

            double Inclination;             // an inclination of the orbit in degrees
            double AscendLong;              // a longitude of the ascending node in degrees
            double Eccentricity;            // an eccentricity of the orbit
            double PeriapsArg;              // a periapsis argument in degrees
            double AvrAnomaly;              // an average anomaly in degrees
            double Frequency;               // a frequency of the rotation (loop per day)
            unsigned int Loop;        // a number of the loop

            TLE() :
                Header{}, Satellite{ 0 }, Classification{ 0 },
                LaunchYear{ 0 }, LaunchNumber{ 0 }, EpochYear{ 0 }, EpochTime{ 0 },
                FirstAcc{ 0 }, SecondAcc{ 0 }, Deceleration{ 0 },
                Inclination{ 0 }, AscendLong{ 0 }, Eccentricity{ 0 },
                PeriapsArg{ 0 }, AvrAnomaly{ 0 }, Frequency{ 0 }, Loop{ 0 }
            {}
            ~TLE() {}

            geometry::Oscul to_oscul(const double mu) const;
        };

        bool load_tle(std::istream& fread, TLE& f);
	}
}