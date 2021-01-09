#pragma once
#include "Atmosphere.h"
#include "GeoPotential.h"
#include "Integrators.h"
#include "Structures.h"
#include "GeneralConstants.h"

namespace ball
{
	using namespace space;
	using namespace geometry;
	using namespace time;
	// Object provides the functionality for ballistic calculation.
	// The trajectory of the center of mass for the space vehicle can be calculated.
	// 
	class Ballistic
	{
	private:
		std::unique_ptr<MultiStepIntegrator<double>> _pMultiStepInt;
		std::unique_ptr<SingleStepIntegrator<double>> _pSingleStepInt;
		std::shared_ptr<IAtmosphere> _pAtmosphere;
		GeoPotential _geoPotential;
		std::vector<std::pair<PV, JD>> _trajectory;
		std::vector<size_t> _loops;
		double _sBall;
		double _deltaTime;
		double _eW, _eFl, _eR;

		PV Function(const PV& vec, const JD& t, const double sb);
		void StartRun(const double startStep, const double continueStep, const size_t index);
		void ContinueRun(const double step, const size_t index);

	public:
		Ballistic(
			const std::shared_ptr<IGravity> pGravity,
			const std::shared_ptr<IAtmosphere> pAtmosphere,
			const size_t harmonics);
		Ballistic(
			const std::shared_ptr<IGravity> pGravity,
			const std::shared_ptr<IAtmosphere> pAtmosphere,
			std::unique_ptr<SingleStepIntegrator<double>> pSingleStepInt,
			std::unique_ptr<MultiStepIntegrator<double>> pMultiStepInt, 
			const size_t harmonics);
		Ballistic(const Ballistic& ball) = delete;
		~Ballistic() {}

		Ballistic& operator = (const Ballistic& ball) = delete;

		// Calculatinf the trajectory.
		// x0 - initial conditions;
		// tk - final time;
		// startStep - period for the single step integrator in seconds;
		// continueStep - period for the multi step integrator in seconds.
		// Throws runtime error when height is out of bounds (MinHeight, MaxHeight).
		void Run(
			const State& x0,
			const JD& tk,
			const double startStep = 5.0,
			const double continueStep = 30.0);
		// Calculating the state parameters.
		// time - time of the trajectory point (should be between the start and final time of the integration);
		// x - will be filled by the content.
		bool GetPoint(const JD& time, State& x) const;
		// Returns the reference to current trajectory.
		const std::vector<std::pair<PV, JD>>& Trajectory() const
		{
			return _trajectory;
		}

		double MinHeight = 1e5, MaxHeight = 1e8;

	};
}
