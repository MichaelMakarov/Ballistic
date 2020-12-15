#pragma once
#include "SingleStepIntegrator.h"
#include "Matrix3x3.h"

namespace ball
{
	class PlahovIntegrator : public SingleStepIntegrator
	{
	private:
		const double _a[3]{ 0.212340538, 0.590533136, 0.911412040 };
		const math::Matrix3x3 _mA;

	public:
		PlahovIntegrator() : _mA(
			math::Matrix3x3::Inv(
				math::Matrix3x3(
					_a[0], _a[0] * _a[0], _a[0] * _a[0] * _a[0],
					_a[1], _a[1] * _a[1], _a[1] * _a[1] * _a[1],
					_a[2], _a[2] * _a[2], _a[2] * _a[2] * _a[2])))
		{}
		~PlahovIntegrator() {}

		types::PV Integrate(const double step) const override;
	};
}