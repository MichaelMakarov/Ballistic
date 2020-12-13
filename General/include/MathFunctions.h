#pragma once
#include "Constants.h"

namespace ball
{
    namespace math
    {
        inline double DegToRad(const double degrees)
        {
            return degrees * PI / 180.0;
        }
        inline double RadToDeg(const double radians)
        {
            return radians * 180.0 / PI;
        }

        long double Factorial(const size_t x);
    }
}