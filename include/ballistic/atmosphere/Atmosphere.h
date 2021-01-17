#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
    namespace space
    {
        using namespace general;
        class IAtmosphere
        {
        public:
            virtual double density(
                const geometry::XYZ& position, 
                const time::JD& time) const = 0;
        };
    }
}