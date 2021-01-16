#pragma once
#include "Geometry.h"
#include "Gravity.h"
#include "Mathematics.h"
#include <memory>

namespace ball
{
	namespace space
	{
		// Class represents a geopotential with its properties
		class GeoPotential
		{
		private:
			std::vector<std::pair<double, double>> _harmonics;
			std::vector<math::LegendreFunction> _functions;
			size_t _count;
			double _eR, _eMu;

		public:
			// Creating the geopotential.
			// pPotential - a pointer to the Earth's gravity model.
			// harmonics - a number of spherical functions in the gravity model to consider.
			GeoPotential(
				const std::shared_ptr<IGravity> pPotential,
				const size_t harmonics);
			GeoPotential(const GeoPotential& gp) = delete;
			~GeoPotential() {}

			GeoPotential& operator = (const GeoPotential& gp) = delete;

			// Calculating the value of Earth's potential
			double operator () (const geometry::RBL& coordinates) const;
			// Calculating the acceleration by the potential
			geometry::XYZ acceleration(const geometry::XYZ& coordinates) const;
		};
	}
}