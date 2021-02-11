#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
    namespace space
    {
        template<class AtmosphereType>
        class IAtmosphere
        {
        public:
            virtual double density(
                const general::math::Vec3& position,
                const general::time::JD& time) const
            {
                return static_cast<const AtmosphereType*>(this)->density(position, time);
            }
        };
    }
}