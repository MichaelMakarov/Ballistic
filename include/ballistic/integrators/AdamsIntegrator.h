#pragma once
#include "Integrators.h"

namespace ball
{
	template<class ... ArgsType>
	class AdamsIntegrator : public MultistepIntegrator<AdamsIntegrator<ArgsType ...>, ArgsType ...>
	{
	private:
		const double _b[8]
		{
			-0.3042245370370370572,
			2.445163690476190421,
			-8.612127976190476986,
			17.37965443121693454,
			-22.02775297619047734,
			18.05453869047619264,
			-9.525206679894179018,
			3.589955357142857295
		};
		const double _c[8]
		{
			0.01136739417989418056,
			-0.09384093915343914849,
			0.343080357142857173,
			-0.732035383597883671,
			1.017964616402116551,
			-1.0069196428571429713,
			1.156159060846560838,
			0.3042245370370370017
		};

	public:
		AdamsIntegrator() : MultistepIntegrator<AdamsIntegrator<ArgsType ...>, ArgsType ...>(8) {}
		~AdamsIntegrator() = default;

		void integrate(
			std::pair<general::math::PV, general::time::JD>* pData,
			const double step,
			general::math::PV& xk,
			general::time::JD& tk,
			const ArgsType& ... args) const
		{
			xk.P1 = xk.P2 = xk.P3 = xk.V1 = xk.V2 = xk.V3 = 0;
			auto xt{ _b[0] * func(pData->first, pData->second, args ...) };
			general::math::PV xb;
			for (size_t i = 1; i < 8; ++i)
			{
				pData++;
				xb = func(pData->first, pData->second, args ...);
				// prediction
				xt += _b[i] * xb;
				// correction
				xk += _c[i - 1] * xb;
			}
			xt *= step;
			xt += pData->first; 
			tk = pData->second;
			tk.add_seconds(static_cast<int>(step));
			xk += _c[7] * func(xt, tk, args ...);
			xk *= step;
			xk += pData->first;
		}
	};
}