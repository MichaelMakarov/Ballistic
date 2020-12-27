#pragma once
#include "PV.h"
#include <functional>

namespace ball
{
	class Integrator
	{
	protected:
		std::function<types::PV(const types::PV&, const double)> _func;

	public:
		Integrator() {}
		~Integrator() {}

		virtual types::PV Integrate(const double step) const = 0;

		void Function(types::PV(* function)(const types::PV&, const double))
		{
			_func = function;
		}
	};
}