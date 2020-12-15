#pragma once
#include "RBL.h"
#include "PotentialModel.h"
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
				const std::unique_ptr<IPotential> pPotential,
				const double eR,
				const double eMu,
				const size_t harmonics);
			~GeoPotential() {}

			unsigned int Harmonics() const { return _count; }

			double operator () (const geometry::RBL& coordinates) const;

		};
	}
}