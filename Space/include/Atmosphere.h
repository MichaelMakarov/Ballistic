#pragma once
#include "Geometry.h"
#include "Times.h"

namespace ball
{
    namespace space
    {
        class IAtmosphere
        {
        public:
            virtual double density(
                const geometry::XYZ& position, 
                const time::JD& time) const = 0;
        };
    }
}