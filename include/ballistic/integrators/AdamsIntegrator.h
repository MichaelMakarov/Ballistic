#pragma once
#include "Integrators.h"

namespace ball
{
	template<Arithmetic R>
	class AdamsIntegrator : public MultistepIntegrator<AdamsIntegrator<R>, R>
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
		AdamsIntegrator() : MultistepIntegrator<AdamsIntegrator<R>, R>(8) {}
		~AdamsIntegrator() = default;

		void integrate(
			std::pair<R, general::time::JD>* pData,
			const double step,
			R& xk,
			general::time::JD& tk) const
		{
			xk.Pos.X = xk.Pos.Y = xk.Pos.Z = xk.Vel.X = xk.Vel.Y = xk.Vel.Z = 0;
			auto xt{ _b[0] * this->func(pData->first, pData->second) };
			general::math::PV xb;
			for (size_t i = 1; i < 8; ++i)
			{
				pData++;
				xb = this->func(pData->first, pData->second);
				// prediction
				xt += _b[i] * xb;
				// correction
				xk += _c[i - 1] * xb;
			}
			xt *= step;
			xt += pData->first; 
			tk = pData->second;
			tk.add_seconds(static_cast<int>(step));
			xk += _c[7] * this->func(xt, tk);
			xk *= step;
			xk += pData->first;
		}
	};
}