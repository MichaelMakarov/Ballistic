#pragma once
#include "Integrator.h"
#include "DateTime.h"
#include <functional>

namespace ball
{
	class SingleStepIntegrator : public Integrator
	{
	protected:
		time::JD _t0;
		types::PV _x0;
		std::function<types::PV(const types::PV&, const time::JD)> _func;

	public:
		using Integrator::Integrator;
		~SingleStepIntegrator() {}

		void Initialize(
			const types::PV& x0,
			const time::JD t0)
		{
			_x0 = x0;
			_t0 = t0;
		}

		void Function(types::PV(* function)(const types::PV&, const time::JD))
		{
			_func = function;
		}
	};
}