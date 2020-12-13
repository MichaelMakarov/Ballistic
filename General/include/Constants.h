#pragma once

namespace ball
{
    namespace time
    {
        inline constexpr unsigned int SecPerMin{ 60 };
        inline constexpr unsigned int SecPerHour{ 3600 };
        inline constexpr unsigned int SecPerDay{ 86400 };

        inline constexpr unsigned int MillisecPerSec{ 1000 };
        inline constexpr unsigned int MillisecPerMin{ 60000 };
        inline constexpr unsigned int MillisecPerHour{ 3600000 };

        inline constexpr unsigned int HoursPerDay{ 24 };

        inline constexpr unsigned int MinPerHour{ 60 };
        inline constexpr unsigned int MinPerDay{ 1440 };

        inline constexpr unsigned int DaysPerYear{ 365 };

        inline constexpr double JD1900{ 2411541.0 };
        inline constexpr double JD2000{ 2451545.0 };

    }

    namespace math
    {
        inline constexpr double PI{ 3.1415926535 };
    }
}