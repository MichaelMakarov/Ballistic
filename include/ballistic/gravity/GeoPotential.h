#pragma once
#include "general/Matrix.h"
#include "EarthModel.h"
#include "general/Mathematics.h"
#include <memory>

namespace ball
{
	// Class represents a geopotential with its properties
	class GeoPotential
	{
	private:
		const std::vector<std::pair<double, double>>& _harmonics;
		std::vector<std::pair<double, double>> _cs;
		std::vector<double> _pnm;
		size_t _count;
		double _eR, _eMu;

		void calc_trigonometric(const double coslambda, const double sinlambda);
		void calc_polynoms(const double cosphi, const double sinphi);
		

	public:
		/// <summary>
		/// Creating the geopotential
		/// </summary>
		/// <param name="eMu"> - a gravitational constant</param>
		/// <param name="eR"> - a radius of equator</param>
		/// <param name="emodel"> - a potential model</param>
		/// <param name="harmonics"> - a number of harmonics</param>
		GeoPotential(
			const double eMu, const double eR,
			const IEarth& emodel,
			const size_t harmonics);
		GeoPotential(const GeoPotential& gp) = delete;
		GeoPotential(GeoPotential&& gp) noexcept;
		~GeoPotential() noexcept = default;

		GeoPotential& operator = (const GeoPotential& gp) = delete;
		GeoPotential& operator = (GeoPotential&& gp) noexcept = default;

		/// <summary>
		/// Calculating the value of Earth's potential
		/// </summary>
		/// <param name="point"> - a point in GCS where to calculate geopotential</param>
		/// <returns>geopotential value</returns>
		double operator () (const general::math::Vec3& point);
		/// <summary>
		/// Calculating the potential acceleration with respect to geocentric cartesian coordinates 
		/// </summary>
		/// <param name="point"> - a point in GCS where to calculate the geopotential</param>
		/// <returns>a vector of partial acceleration (dU/dx, dU/dy, dU/dz)</returns>
		general::math::Vec3 acceleration(const general::math::Vec3& point);
		/// <summary>
		/// Calculating the potential acceleration of the first and second degree with respect to geocentric cartesian coordinates
		/// </summary>
		/// <param name="point"> - a point in GCS where to calculate the geopotential</param>
		/// <returns>a vector of first degree acceleration and matrix of second degree partial acceleration</returns>
		std::pair<general::math::Vec3, general::math::Matrix3x3> derivatives(const general::math::Vec3& point);
	};
}