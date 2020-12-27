#pragma once
#include "Integrator.h"
#include <vector>

namespace ball
{
	class MultiStepIntegrator : public Integrator
	{
	protected:
		std::vector<std::pair<types::PV, double>> _xList;

	private:
		const size_t _degree;

	public:
		MultiStepIntegrator(const size_t degree) : Integrator(), _degree(degree) {}
		~MultiStepIntegrator() {}

		bool Initialize(const std::initializer_list<std::pair<types::PV, double>> xList)
		{
			if (xList.size() == _degree)
			{
				_xList = xList;
				return true;
			}
			return false;
		}

		size_t Degree() const
		{
			return _degree;
		}
	};
}