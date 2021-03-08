#pragma once
#include "general/Geometry.h"
#include "EarthModel.h"
#include "general/Mathematics.h"
#include <memory>

namespace ball
{
	// Class represents a geopotential with its properties
	class GeoPotential
	{
	private:
		std::vector<std::pair<double, double>> _harmonics;
		std::vector<std::pair<double, double>> _cs;
		std::vector<double> _pnm;
		size_t _count;
		double _eR, _eMu;

		void calc_trigonometric(const double coslambda, const double sinlambda);
		void calc_polynoms(const double cosphi, const double sinphi);

	public:
		// Creating the geopotential.
		// pPotential - a pointer to the Earth's gravity model.
		// harmonics - a number of spherical functions in the gravity model to consider.
		GeoPotential(
			const IEarth& emodel,
			const size_t harmonics);
		GeoPotential(const GeoPotential& gp);
		GeoPotential(GeoPotential&& gp) noexcept;
		~GeoPotential() noexcept = default;

		GeoPotential& operator = (const GeoPotential& gp);
		GeoPotential& operator = (GeoPotential&& gp) noexcept;

		// Calculating the value of Earth's potential
		double operator () (const general::math::Vec3& coordinates);
		// Calculating the acceleration by the potential
		general::math::Vec3 acceleration(const general::math::Vec3& coordinates);
	};
}