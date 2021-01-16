#pragma once
#include "Atmosphere.h"
#include <vector>

namespace ball
{
    namespace space
    {
        class StaticAtmosphere81 : public IAtmosphere
        {
        public:
            StaticAtmosphere81(
                const double eR, 
                const double eFl) : _eR{ eR }, _eFl{ eFl } {}
            ~StaticAtmosphere81() {}

            double density(
                const geometry::XYZ& position,
                const time::JD& time) const override;
        private:
            const double _eR, _eFl;

            static const std::vector<double> _height;
            static const std::vector<double> _a0;
            static const std::vector<double> _k1;
            static const std::vector<double> _k2;
        };
    }
}