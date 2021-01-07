#pragma once
#include "PV.h"
#include "DateTime.h"
#include <functional>

namespace ball
{
	template<class ... ArgType>
	class Integrator
	{
	public:
		Integrator() {}
		~Integrator() {}

		virtual types::PV Integrate(const double step, const ArgType ... args) const = 0;

		std::function<types::PV(const types::PV&, const time::JD& t, const ArgType ... args)> Func;
	};
}