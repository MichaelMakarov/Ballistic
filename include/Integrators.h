#pragma once
#include "Geometry.h"
#include "Times.h"
#include <functional>

namespace ball
{
	template<class ... ArgType>
	class Integrator
	{
	public:
		Integrator() {}
		~Integrator() {}

		virtual geometry::PV Integrate(const double step, const ArgType ... args) const = 0;

		std::function<geometry::PV(const geometry::PV&, const time::JD& t, const ArgType ... args)> Func;
	};

	template<class ... ArgType>
	class SingleStepIntegrator : public Integrator<ArgType ...>
	{
	protected:
		time::JD _t0;
		geometry::PV _x0;

	public:
		SingleStepIntegrator() : _t0(0), _x0() {}
		~SingleStepIntegrator() {}

		void Initialize(
			const geometry::PV& x0,
			const time::JD& t0)
		{
			_x0 = x0;
			_t0 = t0;
		}
	};

	template<class ... ArgType>
	class MultiStepIntegrator : public Integrator<ArgType ...>
	{
	protected:
		std::pair<geometry::PV, time::JD>* _pData;

	private:
		const size_t _degree;

	public:
		explicit MultiStepIntegrator(const size_t degree = 1) :	_degree(degree), _pData{ nullptr } {}
		~MultiStepIntegrator() {}

		template<class Iterator>
		void Initialize(Iterator iterData)
		{
			_pData = iterData._Ptr;
		}

		void Initialize(std::pair<geometry::PV, time::JD>* const pData)
		{
			_pData = pData;
		}

		size_t Degree() const
		{
			return _degree;
		}
	};
}