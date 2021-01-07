#pragma once
#include "Integrator.h"

namespace ball
{
	template<class ... ArgType>
	class SingleStepIntegrator : public Integrator<ArgType ...>
	{
	protected:
		time::JD _t0;
		types::PV _x0;

	public:
		SingleStepIntegrator() : _t0(0), _x0() {}
		~SingleStepIntegrator() {}

		void Initialize(
			const types::PV& x0,
			const time::JD& t0)
		{
			_x0 = x0;
			_t0 = t0;
		}
	};
}