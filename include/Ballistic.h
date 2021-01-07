#pragma once
#include "AtmosphereModel.h"
#include "GeoPotential.h"
#include "SingleStepIntegrator.h"
#include "MultiStepIntegrator.h"
#include "StateTypes.h"
#include "Constants.h"

namespace ball
{
	using namespace tasks;
	using namespace types;
	using namespace time;
	class Ballistic
	{
	private:
		std::unique_ptr<MultiStepIntegrator<double>> _pMultiStepInt;
		std::unique_ptr<SingleStepIntegrator<double>> _pSingleStepInt;
		std::unique_ptr<IAtmosphere> _pAtmosphere;
		std::unique_ptr<GeoPotential> _pGeoPotential;
		std::vector<std::pair<PV, JD>> _trajectory;
		std::vector<size_t> _loops;
		double _sBall;
		double _step;
		double _eW, _eFl, _eR;

		PV Function(const PV& vec, const JD& t, const double sb);
		void StartRun(const double startStep, const double continueStep, const size_t index);
		void ContinueRun(const double step, const size_t index);

	public:
		Ballistic() :	
			_step{ 0 },
			_sBall{ 0 },
			_eW{ 0 }, _eFl{ 0 }, _eR{ 0 }{}
		Ballistic(
			std::unique_ptr<IGravity> pGravity,
			std::unique_ptr<IAtmosphere> pAtmosphere,
			const size_t harmonics);
		Ballistic(
			std::unique_ptr<IGravity> pGravity,
			std::unique_ptr<IAtmosphere> pAtmosphere,
			std::unique_ptr<SingleStepIntegrator<double>> pSingleStepInt,
			std::unique_ptr<MultiStepIntegrator<double>> pMultiStepInt, 
			const size_t harmonics);
		Ballistic(const Ballistic& ball) = delete;
		~Ballistic() {}

		Ballistic& operator = (const Ballistic& ball) = delete;

		void Run(
			const StateParams& x0,
			const JD& tk,
			const double startStep = 5.0,
			const double continueStep = 30.0);

		const std::vector<std::pair<PV, JD>>& Trajectory() const
		{
			return _trajectory;
		}


		double MinHeight = 1e5, MaxHeight = 1e8;

	};
}
