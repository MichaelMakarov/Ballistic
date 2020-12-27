#pragma once
#include "MultiStepIntegrator.h"

namespace ball
{
	class AdamsBashforthIntegrator : public MultiStepIntegrator
	{
	private:
		const double _b[8]{ 3.5975, -9.5252, 18.0545, -22.0278, 17.3797, -8.6121, 2.4452, -0.3044 };

	public:
		AdamsBashforthIntegrator() : MultiStepIntegrator(8) {}
		~AdamsBashforthIntegrator() {}

		types::PV Integrate(const double step) const override;
	};
}