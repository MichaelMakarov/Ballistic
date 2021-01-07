#pragma once
#include "Integrator.h"
#include <vector>
#include <list>

namespace ball
{
	template<class ... ArgType>
	class MultiStepIntegrator : public Integrator<ArgType ...>
	{
	protected:
		std::pair<types::PV, time::JD>* _pData;

	private:
		const size_t _degree;

	public:
		explicit MultiStepIntegrator(const size_t degree = 1) : 
			Integrator(), 
			_degree(degree), 
			_pData{ nullptr }
		{}
		~MultiStepIntegrator() {}

		template<class Iterator>
		void Initialize(Iterator iterData)
		{
			_pData = iterData._Ptr;
		}

		void Initialize(const std::vector<std::pair<types::PV, time::JD>>::iterator iter)
		{
			_pData = iter._Ptr;
		}
		void Initialize(const std::initializer_list<std::pair<types::PV, time::JD>>::iterator iter)
		{
			_pData = iter;
		}
		

		void Initialize(std::pair<types::PV, time::JD>* pData)
		{
			_pData = pData;
		}

		size_t Degree() const
		{
			return _degree;
		}
	};
}