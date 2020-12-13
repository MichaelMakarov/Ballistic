#pragma once
#include "StateTypes.h"
#include <string>

namespace ball
{
    namespace formats
    {
        struct TleFormat
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
            unsigned int LoopNumber;        // a number of the loop

            TleFormat() {}
            ~TleFormat() {}

            types::OsculParams ToOsculParams() const;

            friend bool TleFromFile (const char* filepath, TleFormat& f);
        };
    }
}