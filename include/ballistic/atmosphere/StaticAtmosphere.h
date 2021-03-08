#pragma once
#include "Atmosphere.h"
#include <vector>

namespace ball
{
    class StaticAtmosphere81 : public IAtmosphere<StaticAtmosphere81>
    {
    public:
        StaticAtmosphere81(
            const double eR,
            const double eFl) : _eR{ eR }, _eFl{ eFl } {}
        ~StaticAtmosphere81() {}

        double density(
            const general::math::Vec3& position,
            const general::time::JD& time) const;
    private:
        const double _eR, _eFl;

        static const std::vector<double> _height;
        static const std::vector<double> _a0;
        static const std::vector<double> _k1;
        static const std::vector<double> _k2;
    };
}