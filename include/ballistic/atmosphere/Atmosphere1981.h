#pragma once
#include "Atmosphere.h"

namespace ball
{
    class Atmosphere1981 : public IAtmosphere<Atmosphere1981>
    {
    public:
        Atmosphere1981(
            const double eR,
            const double eFl) : _eR{ eR }, _eFl{ eFl } {}
        ~Atmosphere1981() = default;

        double density(
            const general::math::Vec3& position,
            const general::time::JD& time) const;
    private:
        const double _eR, _eFl;

        static const double _height[9];
        static const double _a0[9];
        static const double _k1[9];
        static const double _k2[9];
    };
}