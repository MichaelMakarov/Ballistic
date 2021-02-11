#pragma once
#include "general/Geometry.h"
#include "EarthModel.h"
#include "general/Mathematics.h"
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
			size_t _count;
			double _eR, _eMu;

		public:
			// Creating the geopotential.
			// pPotential - a pointer to the Earth's gravity model.
			// harmonics - a number of spherical functions in the gravity model to consider.
			GeoPotential(
				const std::shared_ptr<IEarth> pPotential,
				const size_t harmonics);
			GeoPotential(const GeoPotential& gp);
			GeoPotential(GeoPotential&& gp) noexcept;
			~GeoPotential() noexcept {}

			GeoPotential& operator = (const GeoPotential& gp);
			GeoPotential& operator = (GeoPotential&& gp) noexcept;

			// Calculating the value of Earth's potential
			double operator () (const general::math::Vec3& coordinates) const;
			// Calculating the acceleration by the potential
			general::math::Vec3 acceleration(const general::math::Vec3& coordinates) const;
		};
	}
}