#pragma once
#include "RBL.h"
#include "XYZ.h"
#include "GravityModel.h"
#include "LegendrePolynom.h"
#include <memory>

namespace ball
{
	namespace tasks
	{
		// an object represents a geopotential with its properties
		class GeoPotential
		{
		private:
			std::vector<std::pair<double, double>> _harmonics;
			std::vector<math::LegendreFunction> _functions;
			size_t _count;
			double _eR, _eMu;



		public:
			GeoPotential(
				const std::unique_ptr<IGravity> pPotential,
				const size_t harmonics);
			GeoPotential(const GeoPotential& gp) = delete;
			~GeoPotential() {}

			GeoPotential& operator = (const GeoPotential& gp) = delete;
			GeoPotential& operator = (GeoPotential&& gp) noexcept;

			size_t Harmonics() const { return _count; }

			double Radius() const { return _eR; }
			double Mu() const { return _eMu; }

			double operator () (const geometry::RBL& coordinates) const;

			geometry::XYZ Acceleration(const geometry::XYZ& coordinates) const;

		};
	}
}