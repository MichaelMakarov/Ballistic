#include "Ballistic.h"
#include "AdamsIntegrator.h"
#include "RungeKuttaIntegrator.h"
#include "AstroValues.h"

namespace ball
{
	PV Ballistic::Function(const PV& vec, const JD& t, const double sb)
	{
		PV func(vec.V1, vec.V2, vec.V3, 0, 0, 0);
		geometry::XYZ xyzPos(vec.P1, vec.P2, vec.P3);
		double r = xyzPos.Length(),
			r_2 = r * r,
			w_2 = _eW * _eW;
		double h = GCS_HeightFromPosition(xyzPos, _eR, _eFl);
		if (h < MinHeight || h > MaxHeight)
			throw std::runtime_error("Height is out of bounds!");
		// geopotential aceleration with a centrifugal and a coriolis force
		auto xyzAcPot{ _pGeoPotential->Acceleration(xyzPos) };
		// atmosphere aceleration a = v * s * rho, 
		// s - a ballistic coefficient,
		// v - a velocity of the vehicle,
		// rho - a density of the atmosphere
		double density = _pAtmosphere->Density(xyzPos, t);
		double v = std::sqrt(vec.V1 * vec.V1 + vec.V2 * vec.V2 + vec.V3 * vec.V3);
		double acAtm = v * density * sb;
		// the addition all the components
		func.V1 = xyzAcPot.X + w_2 * vec.P1 + 2 * _eW * vec.V2 - acAtm * vec.V1;
		func.V2 = xyzAcPot.Y + w_2 * vec.P2 - 2 * _eW * vec.V1 - acAtm * vec.V2;
		func.V3 = xyzAcPot.Z - acAtm * vec.V3;
		return func;
	}

	void Ballistic::StartRun(
		const double startStep, 
		const double continueStep, 
		const size_t index)
	{
		std::pair<PV, JD> x;
		PV vec;
		size_t loop;
		bool intersection;
		const size_t n{ static_cast<size_t>(continueStep / startStep) };
		const double step{ continueStep / n };
		const double dt{ step / time::SEC_PER_DAY };
		for (size_t i = 0; i < index; ++i)
		{
			x = _trajectory[i];
			loop = _loops[i];
			for (size_t k = 0; k < n; ++k)
			{
				_pSingleStepInt->Initialize(x.first, x.second);
				vec = _pSingleStepInt->Integrate(step, _sBall);
				intersection = std::signbit(x.first.P3) == false && std::signbit(vec.P3) == true;
				x = { vec, x.second + dt };
				loop = loop - intersection ? 1 : 0;
			}
			_trajectory[i + 1] = x;
			_loops[i + 1] = loop;
		}
	}

	void Ballistic::ContinueRun(const double step, const size_t index)
	{
		std::pair<PV, JD> x;
		PV vec;
		bool intersection;
		const double dt{ step / time::SEC_PER_DAY };
		auto iter{ _trajectory.begin() };
		for (size_t i = index; i < _trajectory.size() - 1; ++i)
		{
			x = _trajectory[i];
			_pMultiStepInt->Initialize(iter);
			vec = _pMultiStepInt->Integrate(step, _sBall);
			intersection = std::signbit(x.first.P3) == true && std::signbit(vec.P3) == false;
			_trajectory[i + 1] = { vec, x.second + dt };
			_loops[i + 1] = _loops[i] + intersection ? 1 : 0;
			iter++;
		}
	}

	Ballistic::Ballistic(
		std::unique_ptr<IGravity> pGravity,
		std::unique_ptr<IAtmosphere> pAtmosphere,
		const size_t harmonics)
	{
		_eR = pGravity->R();
		_eFl = pGravity->Fl();
		_eW = pGravity->W();
		_pGeoPotential = { std::make_unique<GeoPotential>(std::move(pGravity), harmonics) };
		_pAtmosphere = std::move(pAtmosphere);
		_pMultiStepInt = { std::make_unique<AdamsIntegrator<double>>() };
		_pSingleStepInt = { std::make_unique<RungeKuttaIntegrator<double>>() };
		auto pFunc = &Ballistic::Function;
		auto pThis = this;
		auto func = [pFunc, pThis](const PV& vec, const JD& time, const double sb) {
			return (pThis->*pFunc)(vec, time, sb);
		};
		_pMultiStepInt->Func = func;
		_pSingleStepInt->Func = func;
	}

	Ballistic::Ballistic(
		std::unique_ptr<IGravity> pGravity,
		std::unique_ptr<IAtmosphere> pAtmosphere,
		std::unique_ptr<SingleStepIntegrator<double>> pSingleStepInt,
		std::unique_ptr<MultiStepIntegrator<double>> pMultiStepInt,
		const size_t harmonics)
	{
		// Earth's constants
		_eR = pGravity->R();
		_eFl = pGravity->Fl();
		_eW = pGravity->W();
		// initializing the physical objects 
		_pGeoPotential = { std::make_unique<GeoPotential>(std::move(pGravity), harmonics) };
		_pAtmosphere = std::move(pAtmosphere);
		// initializing the integrators
		_pSingleStepInt = std::move(pSingleStepInt);
		_pMultiStepInt = std::move(pMultiStepInt);
		auto pFunc = &Ballistic::Function;
		auto pThis = this;
		auto func = [pFunc, pThis](const PV& vec, const JD& time, const double sb) {
			return (pThis->*pFunc)(vec, time, sb);
		};
		_pMultiStepInt->Func = func;
		_pSingleStepInt->Func = func;
	}

	void Ballistic::Run(
		const StateParams& x0,
		const time::JD& tk,
		const double startStep,
		const double continueStep)
	{
		auto second{ 1.0 / time::SEC_PER_DAY };
		if (startStep < second || continueStep < second)
			throw std::invalid_argument("Invalid time step (should be > 1.0 sec)!");
		if (startStep > continueStep)
			throw std::invalid_argument("Invalid start step > continue step!");
		if (tk <= x0.T)
			throw std::invalid_argument("Invalid tk < tn!");
		_step = continueStep;
		_sBall = x0.Sb;
		size_t index = _pMultiStepInt->Degree() - 1;
		size_t count = static_cast<size_t>((tk - x0.T) / _step * time::SEC_PER_DAY) + 1;
		_trajectory.resize(count);
		_loops.resize(count);

		_trajectory[0] = { x0.Vec, x0.T };
		_loops[0] = x0.LoopN;
		// performing initial calculations
		StartRun(startStep, continueStep, index);
		// performing remaining calculations
		ContinueRun(continueStep, index);
	}
}

