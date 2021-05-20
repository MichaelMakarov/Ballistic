#pragma once
#include "Integrators.h"

namespace ball
{
	template<Arithmetic R, Time T>
	class AdamsIntegrator : public MultistepIntegrator<AdamsIntegrator<R, T>, R, T>
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
		AdamsIntegrator() : MultistepIntegrator<AdamsIntegrator<R, T>, R, T>(8) {}
		~AdamsIntegrator() = default;

		template<class Inv>
		void integrate(
			const R* pR, const T* pT,
			const double step,
			R& xk,
			T& tk,
			const Func<Inv, R, T>& func) const
		{
			xk = R();
			auto xt{ func(*pR, *pT) * _b[0] };
			R xb;
			for (size_t i = 1; i < 8; ++i) {
				++pR;
				++pT;
				xb = func(*pR, *pT);
				xt += xb * _b[i];		// prediction
				xk += xb * _c[i - 1];	// correction
			}
			xt *= step;
			xt += (*pR); 
			tk = (*pT) + step;
			xk += func(xt, tk) * _c[7];
			xk *= step;
			xk += (*pR);
		}
	};
}