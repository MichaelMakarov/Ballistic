#include "GeoPotential.h"
#include <cmath>

namespace ball
{
	namespace tasks
	{
		GeoPotential::GeoPotential(
			const std::unique_ptr<IPotential> pPotential,
			const double eR,
			const double eMu,
			const size_t count)
		{
			_count = count > pPotential->Count() ? pPotential->Count() : count;
			_eR = eR;
			_eMu = eMu;
			size_t dim = 0;
			for (size_t i = 0; i <= _count; ++i) dim += i + 1;
			auto harmonics = pPotential->Harmonics();
			_functions.resize(dim);
			_harmonics.resize(dim);
			dim = 0;
			for (size_t i = 0; i <= _count; ++i)
				for (size_t j = 0; j <= i; ++j)
				{
					_harmonics[dim] = harmonics[dim];
					_functions[dim++] = math::LegendreFunction(i, j);
				}
		}

		double GeoPotential::operator () (const geometry::RBL& coordinates) const
		{
			double result = 0.0,
				sinphi = std::sin(coordinates.B),
				R_r = _eR / coordinates.R,
				l = coordinates.L,
				mult = 1,
				p;
			size_t k = 0;
			std::vector<std::pair<double, double>> ls(_count + 1);
			for (size_t i = 0; i <= _count; ++i) 
				ls[i] = std::pair<double, double>(std::cos(i * l), std::sin(i * l));
			for (size_t n = 0; n <= _count; ++n)
			{
				for (size_t m = 0; m <= n; ++m)
				{
					p = _harmonics[k].first * ls[m].first + _harmonics[k].second * ls[m].second;
					result += mult * p * _functions[k](sinphi) * _functions[k].Norm();
					k++;
				}
				mult *= R_r;
			}
			return _eMu / coordinates.R * result;
		}
	}
}