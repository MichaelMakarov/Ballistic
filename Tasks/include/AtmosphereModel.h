#pragma once
#include "XYZ.h"
#include "DateTime.h"

namespace ball
{
    namespace tasks
    {
        class IAtmosphere
        {
        public:
            virtual double Density(
                const geometry::XYZ& position, 
                const time::JD& time) = 0;
        };
    }
}