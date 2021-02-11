#pragma once
#include "Forecast.h"
#include "gravity/GeoPotential.h"
#include "atmosphere/Atmosphere.h"
#include "Conversions.h"

namespace ball
{
	// A forecast implements atmosphere standard and geopotential model
	template<class AtmType>
	class PAForecast : public Forecast<PAForecast<AtmType>>
	{
	private:
		GeoPotential _geopotential;
		std::shared_ptr<IAtmosphere<AtmType>> _pAtmosphere;
		double _eW, _eFl, _eR;

	public:
		PAForecast(
			const std::shared_ptr<IEarth>& pGravity,
			const size_t harmonics,
			const std::shared_ptr<IAtmosphere<AtmType>>& pAtmosphere) :
			_pAtmosphere{ pAtmosphere },
			_geopotential{ GeoPotential(pGravity, harmonics) },
			_eW{ pGravity->W() },
			_eFl{ pGravity->Fl() },
			_eR{ pGravity->R() }
		{}
		PAForecast(const PAForecast& f) = default;
		PAForecast(PAForecast&& f) noexcept : Forecast(f),
			_geopotential{ std::move(f._geopotential) },
			_pAtmosphere{ std::move(f._pAtmosphere) },
			_eW{ f._eW }, _eFl{ f._eFl }, _eR{ f._eR }
		{}
		~PAForecast() = default;

		PAForecast& operator = (const PAForecast& f) = default;
		PAForecast& operator = (PAForecast&& f) noexcept
		{
			_geopotential = std::move(f._geopotential);
			_pAtmosphere = std::move(f._pAtmosphere);
			_eW = f._eW;
			_eFl = f._eFl;
			_eR = f._eR;
			_eW = _eR = _eFl = 0;
			return *this;
		}

		// Acelerations calculation using current vector in GCS and time
		general::math::PV function(const general::math::PV& vec, const general::time::JD& t)
		{
			general::math::PV ac(vec.V1, vec.V2, vec.V3, 0, 0, 0);
			auto xyzPos{ general::math::Vec3(vec.P1, vec.P2, vec.P3) };
			double r = xyzPos.length();
			double r_2{ r * r };
			double w_2{ _eW * _eW };
			double h = GCS_height_from_position(xyzPos, _eR, _eFl);
			if (h < MinHeight || h > MaxHeight)
				throw std::runtime_error("Height is out of bounds!");
			// geopotential aceleration with a centrifugal and a coriolis force
			auto xyzAcPot{ _geopotential.acceleration(xyzPos) };
			// atmosphere aceleration a = v * s * rho, 
			// s - a ballistic coefficient,
			// v - a velocity of the vehicle,
			// rho - a density of the atmosphere
			double density = _pAtmosphere->density(xyzPos, t);
			double v = std::sqrt(vec.V1 * vec.V1 + vec.V2 * vec.V2 + vec.V3 * vec.V3);
			double acAtm = v * density * sBall;
			// the addition all the components
			ac.V1 = xyzAcPot.X + w_2 * vec.P1 + 2 * _eW * vec.V2 - acAtm * vec.V1;
			ac.V2 = xyzAcPot.Y + w_2 * vec.P2 - 2 * _eW * vec.V1 - acAtm * vec.V2;
			ac.V3 = xyzAcPot.Z - acAtm * vec.V3;
			return ac;
		}
	};
}