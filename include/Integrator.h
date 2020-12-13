#pragma once
#include "PV.h"
#include "DateTime.h"


namespace ball
{
	class Integrator
	{
	public:
		Integrator() {}
		~Integrator() {}

		virtual types::PV Integrate(const double step) const = 0;
	};
}