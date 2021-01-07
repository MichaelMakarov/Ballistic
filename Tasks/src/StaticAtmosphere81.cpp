#include "StaticAtmosphere81.h"
#include "AstroValues.h"

namespace ball
{
	namespace tasks
	{
		double StaticAtmosphere81::Density(
			const geometry::XYZ& position,
			const time::JD& time) const
		{
			double h = GCS_HeightFromPosition(position, _eR, _eFl);
			size_t i;
			for (i = 0; i < 9; ++i)
				if (h < _height[i]) break;
			h -= _height[--i];
			return _a0[i] * std::exp(h * (h * _k1[i] - _k2[i]));
		}
	}
}