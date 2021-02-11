#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
    namespace space
    {
        using namespace general;
        template<class AtmosphereType>
        class IAtmosphere
        {
        public:
            virtual double density(
                const geometry::XYZ& position,
                const time::JD& time) const
            {
                return static_cast<const AtmosphereType*>(this)->density(position, time);
            }
        };
    }
}