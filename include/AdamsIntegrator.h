#pragma once
#include "MultiStepIntegrator.h"

namespace ball
{
	template<class ... ArgType>
	class AdamsIntegrator : public MultiStepIntegrator<ArgType ...>
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
		AdamsIntegrator() : MultiStepIntegrator(8) {}
		~AdamsIntegrator() {}

		types::PV Integrate(const double step, const ArgType ... args) const override
		{
			types::PV x;
			time::JD t;
			std::pair<PV, JD>* pData = _pData;
			// prediction
			for (size_t i = 0; i < 8; ++i)
			{
				x += _b[i] * Func(pData->first, pData->second, args ...);
				pData++;
			}
			x *= step;
			pData--;
			x += pData->first;
			// correction
			t = pData->second;
			t.AddSeconds(static_cast<int>(step));
			x = _c[7] * Func(x, t, args ...);
			pData -= 6;
			for (size_t i = 0; i < 7; ++i)
			{
				x += _c[i] * Func(pData->first, pData->second, args ...);
				pData++;
			}
			x *= step;
			pData--;
			x += pData->first;
			return x;
		}
	};
}