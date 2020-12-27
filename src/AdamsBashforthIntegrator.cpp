#include "AdamsBashforthIntegrator.h"

namespace ball
{
	types::PV AdamsBashforthIntegrator::Integrate(const double h) const
	{
		types::PV result;
		for (size_t i = 0; i < 8; ++i)
		{
			result += _b[i] * _func(_xList[i].first, _xList[i].second);
		}
		result *= h;
		result += _xList[7].first;
		return result;
	}
}