#pragma once
#include "Integrator.h"
#include <functional>

namespace ball
{
	class SingleStepIntegrator : public Integrator
	{
	protected:
		double _t0;
		types::PV _x0;

	public:
		using Integrator::Integrator;
		~SingleStepIntegrator() {}

		void Initialize(
			const types::PV& x0,
			const double t0)
		{
			_x0 = x0;
			_t0 = t0;
		}
	};
}