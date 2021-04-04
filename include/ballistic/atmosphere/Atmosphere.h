#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
    template<class Atmosphere, class ... Args>
    class IAtmosphere
    {
    public:
        double density(
            const general::math::Vec3& position,
            const general::time::JD& time,
            const Args& ... args) const
        {
            return static_cast<const Atmosphere*>(this)->density(position, time, args ...);
        }
    };
}